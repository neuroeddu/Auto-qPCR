import re


def rx_rename(col_name):
	# compiles regular expressions
	rct = re.compile('((?<![\w _])[(]*c[()ycle ]*t[)hreshold]*(?![\w\W]))', re.IGNORECASE)
	rquant = re.compile('((?<![\w _])[(]*quant[ity)]*\Z(?! sd)(?! mean))|(?<![\w _])(ng)\Z(?! sd)', re.IGNORECASE)
	rsamp = re.compile('(samp)+|(chrom)+|(desc)+', re.IGNORECASE)
	rtarg = re.compile('(targ)+|(gene)+', re.IGNORECASE)
	rdye = re.compile('(repor)+|(dye)+|(fluor)+', re.IGNORECASE)
	rtask = re.compile('(task)+|(role)+|(content)+', re.IGNORECASE)

	if re.match(rct, col_name):
		return str('CT')
	if re.match(rquant, col_name):
		return str('Quantity')
	if re.match(rsamp, col_name):
		return str('Sample Name')
	if re.match(rtarg, col_name):
		return str('Target Name')
	if re.match(rdye, col_name):
		return str('Reporter')
	if re.match(rtask, col_name):
		return str('Task')
	else:
		return col_name