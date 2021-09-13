# Auto-q-PCR is a program for analysis of qPCR data for absolute and relative quantification
# Copyright (C) 2021 Rhalena Thomas, Eddie Cai, Gracia Gu and Iva Demirova
#This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation,  version 3 of the License.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.


#########################################################################################################

# called to make inputs case insensitive

import re


def rx_rename(col_name):
	# compiles regular expressions
	rct = re.compile('((?<![\w _])[(]*c[()ycle ]*t[)hreshold]*(?![\w\W]))|CQ', re.IGNORECASE)
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
