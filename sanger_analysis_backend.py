import sys
import os
from io import BytesIO
import matplotlib.pylab as plt
import base64


sys.path.append('%s/%s' % (os.path.dirname(os.path.realpath(__file__)), 'BIU'))
import biu

hg = biu.genomes.Ensembl(grch37=True)

def fig_to_uri(in_fig, close_all=True, **save_args):
    # type: (plt.Figure) -> str
    """
    Save a figure as a URI
    :param in_fig:
    :return:
    """
    out_img = BytesIO()
    biu.utils.fs.mkdirp('images')
    in_fig.savefig('images/image.png', format='png', **save_args)
    in_fig.savefig(out_img, format='png', **save_args)
    if close_all:
        in_fig.clf()
        plt.close('all')
    out_img.seek(0)  # rewind file
    encoded = base64.b64encode(out_img.read()).decode("ascii").replace("\n", "")
    return "data:image/png;base64,{}".format(encoded)
#edef


# Determine if the primer matches to the correct variantID
def _var_is_in_region(sseqid, sstart, send, variantid):
	"""
	Determine if a region defined by sseqid:sstart-send contains variantid
	Inputs:
	  sseqid : String. chromosomeID
	  sstart : Integer. Start of region
	  send : Integer. End of region
	  variantID : String. Format ChromosomeID-Position-Reference-Alternative

	Output:
	  Is variant in region? Boolean.
	"""
	var_chr, var_pos, var_ref, var_alt = variantid.split('-')
	var_pos = int(var_pos)
	sstart, send = min(sstart, send), max(sstart, send)
	if str(sseqid).isdigit():
		sseqid = '%02d' % int(sseqid)
	#fi
	if (sseqid == var_chr) and (var_pos > int(sstart)) and (var_pos < int(send)):
		return True
	#fi
	return False
#edef

def _primer_true_vid(primers_table, sseqid, sstart, send):

	match = []
	for (i, (primer, variantno, variantid)) in primers_table[["primer", "variantno", "variantid"]].iterrows():
		if _var_is_in_region(sseqid, sstart, send, variantid):
			match.append(variantno)
		#fi
	return set(match)
#edef

################################################################################################################

def determine_primer_alignments(primers):
	"""
	determine primer alignments
	Input:
	  primers: DataFrame with columns: ...
	Output:
	  bres : Dataframe with results from Blast and true alignments
	  incorrect: Dict with incorrect primer mappings
	"""
	primer_fasta = [ biu.formats.Sequence(str(p['primer']) + '-' + str(p['direction']), str(p['sequence']))
				 for (i,p) in primers.iterrows() ]
	primer_fasta = biu.formats.Fasta(primer_fasta)

	bres = biu.pipelines.Blast(hg.genome, primer_fasta,
						   config={"e_threshold" : 10,
								   "options":"-task blastn-short"}# We must specify that we have short sequences
						  ).getResult()

	# Select the hit with the best match
	bres = bres[bres.evalue == bres.groupby('qseqid').transform(min).evalue]
	# Split primer and direction
	bres['primer'], bres['dir'] = zip(*bres.qseqid.apply(lambda s: s.split('-')))
	# Groupby primer, and get the full range of where the primer pairs should align
	bres = bres.groupby('primer').agg(list)[["sseqid", "sstart", "send"]]
	bres['sseqid'] = bres['sseqid'].apply(lambda s: s[0])
	bres['sstart'] = bres['sstart'].apply(lambda s: min(s))
	bres['send'] = bres['send'].apply(lambda s: max(s))
	# Join in the original variant information
	bres = bres.join(primers.set_index('primer')[['variantno', 'variantid']])

	bres['correct'] = bres.apply(lambda r: _var_is_in_region(*list(r[["sseqid", "sstart", "send", 'variantid']])),
								 axis=1)

	bres['true_vid'] = bres.apply(lambda r: _primer_true_vid(primers, *list(r[["sseqid", "sstart", "send"]])), axis=1)

	incorrect = {}
	for i, row in bres[~bres.correct].drop_duplicates('variantno').iterrows():
		incorrect[row.variantid] = row.true_vid
		print('%3s -> ' % row.variantno, row.true_vid)
	#efor
	return bres, incorrect
#edef

def align_abi_files(abi_dict):
	config = { "blast_fields" : "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen sseq qseq" }
	traces = biu.formats.Fasta([ biu.formats.Sequence(n, abi_dict[n].seq) for n in abi_dict])
	bres_abi = biu.pipelines.Blast(hg.genome, traces, config=config).getResult()

	return bres_abi
#edef

def draw_chromatogram(abi, **kwargs):
	ax = abi.chromatogram(**kwargs)
	return ax
#edef

def _determine_where_mutation_should_be(primers, qstart, qend, qseq, sstart, send, sseq, true_vid, well):
    chrpos = {}
    #sstart, send = min(sstart, send), max(sstart, send)
    for tv in true_vid:
        chr, pos, ref, alt = primers[primers.variantno == tv].variantid.values[0].split('-')
        pos = int(pos)

        if sstart < send:
            first_try = qstart + (pos - sstart)
            ngaps_to_pos_in_sseq = len([ c for c in sseq[:first_try] if c == '-' ])
            ngaps_to_pos_in_qseq = len([ c for c in qseq[:first_try] if c == '-' ])
            chrpos[tv] = first_try - ngaps_to_pos_in_qseq + ngaps_to_pos_in_sseq
        else:
            first_try = qend - (pos-send)
            ngaps_to_pos_in_sseq = len([ c for c in sseq[-first_try:] if c == '-' ])
            ngaps_to_pos_in_qseq = len([ c for c in qseq[-first_try:] if c == '-' ])
            chrpos[tv] = first_try - ngaps_to_pos_in_qseq + ngaps_to_pos_in_sseq
    #efor
    return chrpos
#edef

def process_abi_alignments(primer_table, well_table, bres_abi):
	bres_abi['well'] = bres_abi.qseqid.apply(lambda x: x.split('_')[0])
	bres_abi['run'] = bres_abi.qseqid.apply(lambda x: x.split('-')[1][:-4])
	bres_abi['id'] = bres_abi.well + '.' + bres_abi.run.map(str)
	bres_abi = bres_abi.join(well_table.set_index('id'), on='id', lsuffix='_orig')
	bres_abi['correct'] = bres_abi.apply(lambda r:
											_var_is_in_region(*list(r[["sseqid", "sstart",
																	  "send", 'variant_id']])),
											axis=1)

	bres_abi['true_vid'] = bres_abi.apply(lambda r: _primer_true_vid(primer_table, *list(r[["sseqid", "sstart", "send"]])), axis=1)

	bres_abi['chromat_pos'] = bres_abi.apply(lambda r:
                           _determine_where_mutation_should_be(primer_table, *list(r[['qstart', 'qend', 'qseq',
                                                                       'sstart', 'send', 'sseq',
                                                                       'true_vid', 'well'
                                                                      ]])),
                                              axis=1)
	return bres_abi
#edef
