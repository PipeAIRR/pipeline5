$HOSTNAME = ""
params.outdir = 'results'  


if (!params.primerSet){params.primerSet = ""} 
if (!params.samples_file){params.samples_file = ""} 
if (!params.mate){params.mate = ""} 
if (!params.reads){params.reads = ""} 
if (!params.smartNNNA_1){params.smartNNNA_1 = ""} 
if (!params.pairawk_bc_script){params.pairawk_bc_script = ""} 
if (!params.smartNNNA_2){params.smartNNNA_2 = ""} 
// Stage empty file to be used as an optional input where required
ch_empty_file_1 = file("$baseDir/.emptyfiles/NO_FILE_1", hidden:true)
ch_empty_file_2 = file("$baseDir/.emptyfiles/NO_FILE_2", hidden:true)

g_3_primerSets_g0_0 = file(params.primerSet, type: 'any')
g_4_sample_file_g0_0 = file(params.samples_file, type: 'any')
Channel.value(params.mate).into{g_6_mate_g_28;g_6_mate_g1_0;g_6_mate_g1_7;g_6_mate_g13_0;g_6_mate_g13_7;g_6_mate_g35_5;g_6_mate_g35_9;g_6_mate_g38_0;g_6_mate_g38_7;g_6_mate_g39_0;g_6_mate_g39_7}
if (params.reads){
Channel
	.fromFilePairs( params.reads , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.set{g_7_reads_g38_0}
 } else {  
	g_7_reads_g38_0 = Channel.empty()
 }

g_25_vprimers_g35_5 = params.smartNNNA_1 && file(params.smartNNNA_1, type: 'any').exists() ? file(params.smartNNNA_1, type: 'any') : ch_empty_file_2
g_29_PAIRAWK_script_g_28 = file(params.pairawk_bc_script, type: 'any')
g_30_cprimers_g35_5 = params.smartNNNA_2 && file(params.smartNNNA_2, type: 'any').exists() ? file(params.smartNNNA_2, type: 'any') : ch_empty_file_1


process generate_barcode_file_generate_barcode_file {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /Barcode_M1S_Z.fasta$/) "barcode_file/$filename"}
input:
 file primer_set from g_3_primerSets_g0_0
 file sample_file from g_4_sample_file_g0_0

output:
 file "Barcode_M1S_Z.fasta"  into g0_0_fasta_file0_g13_0, g0_0_fasta_file0_g1_0, g0_0_fasta_file0_g38_0, g0_0_fasta_file0_g39_0

shell:
	
'''
awk -F '\t' 'NR==FNR && NR>1 {M1S[$1]=$2; Z[$1]=$3;} NR!=FNR  {print (">"$1"-M1S"); print (M1S[$2]); print (">"$1"-Z"); print (Z[$2]);} ' !{primer_set} !{sample_file} > Barcode_M1S_Z.fasta
'''
}


process mask_primer_m1s_z_start_4_Mask_M1S_and_Z_primers {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_primers-pass.fastq$/) "m1s_z_start_4_pass/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /MP_.*$/) "MP_log_start_4/$filename"}
input:
 set val(name),  file(reads) from g_7_reads_g38_0
 val mate from g_6_mate_g38_0
 file primers from g0_0_fasta_file0_g38_0

output:
 set val(name), file("*_primers-pass.fastq")  into g38_0_reads0_g_23
 file "MP_*"  into g38_0_logFile1_g38_7
 set val(name), file("*_primers-fail.fastq") optional true  into g38_0_reads_failed2_g39_0

script:
method = params.mask_primer_m1s_z_start_4_Mask_M1S_and_Z_primers.method
maxerror = params.mask_primer_m1s_z_start_4_Mask_M1S_and_Z_primers.maxerror
revpr = params.mask_primer_m1s_z_start_4_Mask_M1S_and_Z_primers.revpr
failed = params.mask_primer_m1s_z_start_4_Mask_M1S_and_Z_primers.failed
pf = params.mask_primer_m1s_z_start_4_Mask_M1S_and_Z_primers.pf
nproc = params.mask_primer_m1s_z_start_4_Mask_M1S_and_Z_primers.nproc
primer_start = params.mask_primer_m1s_z_start_4_Mask_M1S_and_Z_primers.primer_start
maxlen = params.mask_primer_m1s_z_start_4_Mask_M1S_and_Z_primers.maxlen
skiprc = params.mask_primer_m1s_z_start_4_Mask_M1S_and_Z_primers.skiprc
//* @style @condition:{method="score",primer_start}, {method="align",maxlen,skiprc} @multicolumn:{method,maxerror,revpr,failed,pf,nrpoc}, {maxlen,skiprc} 

readArray = reads.toString().split(' ')	
revpr_arg = (revpr=="true") ? "--revpr" : ""
failed_arg = (failed=="true") ? "--failed" : "" 

if(method=="score"){
	start = "--start ${primer_start}"
}else{
	start = ""
}

if(method=="align"){
	skiprc_arg = (skiprc=="true") ? "--skiprc" : ""
	maxlen_arg = "--maxlen ${maxlen}"
}else{
	skiprc_arg = ""
	maxlen_arg = ""
}

pf = "--pf ${pf}"

if(mate=="pair"){
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	R1_name =  R1 - '.fastq' - '.fastq'
	R2_name =  R2 - '.fastq' - '.fastq'
	"""
	MaskPrimers.py ${method} -s $R1 -p ${primers} ${start} --mode cut --maxerror ${maxerror} --outname ${R1_name} --log MP_${R1_name}.log  --nproc ${nproc} ${revpr_arg} ${failed_arg} ${skiprc_arg} ${maxlen_arg} ${pf}
	MaskPrimers.py ${method} -s $R2 -p ${primers} ${start} --mode cut --maxerror ${maxerror} --outname ${R2_name} --log MP_${R2_name}.log  --nproc ${nproc} ${revpr_arg} ${failed_arg} ${skiprc_arg} ${maxlen_arg} ${pf}
	"""
}else{
	R1 = readArray[0]
	"""
	echo -e "Assuming inputs for R1"
	MaskPrimers.py ${method} -s $R1 -p ${primers} ${start} --mode cut --maxerror ${maxerror} --outname ${name}_R1 --log MP_R1_${name}.log  --nproc ${nproc} ${revpr_arg} ${failed_arg} ${skiprc_arg} ${maxlen_arg} ${pf}
	"""
}



}


process mask_primer_m1s_z_start_4_parse_log_MP {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tab$/) "m1s_z_log_start_4_tab/$filename"}
input:
 val mate from g_6_mate_g38_7
 file log_file from g38_0_logFile1_g38_7

output:
 file "*.tab"  into g38_7_logFile00

script:
readArray = log_file.toString()	

"""
ParseLog.py -l ${readArray}  -f ID PRIMER BARCODE ERROR
"""

}


process mask_primer_m1s_z_start_5_Mask_M1S_and_Z_primers {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_primers-pass.fastq$/) "m1s_z_start_5_pass/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /MP_.*$/) "MP_log_start_5/$filename"}
input:
 set val(name),  file(reads) from g38_0_reads_failed2_g39_0
 val mate from g_6_mate_g39_0
 file primers from g0_0_fasta_file0_g39_0

output:
 set val(name), file("*_primers-pass.fastq")  into g39_0_reads0_g_23
 file "MP_*"  into g39_0_logFile1_g39_7
 set val(name), file("*_primers-fail.fastq") optional true  into g39_0_reads_failed2_g1_0

script:
method = params.mask_primer_m1s_z_start_5_Mask_M1S_and_Z_primers.method
maxerror = params.mask_primer_m1s_z_start_5_Mask_M1S_and_Z_primers.maxerror
revpr = params.mask_primer_m1s_z_start_5_Mask_M1S_and_Z_primers.revpr
failed = params.mask_primer_m1s_z_start_5_Mask_M1S_and_Z_primers.failed
pf = params.mask_primer_m1s_z_start_5_Mask_M1S_and_Z_primers.pf
nproc = params.mask_primer_m1s_z_start_5_Mask_M1S_and_Z_primers.nproc
primer_start = params.mask_primer_m1s_z_start_5_Mask_M1S_and_Z_primers.primer_start
maxlen = params.mask_primer_m1s_z_start_5_Mask_M1S_and_Z_primers.maxlen
skiprc = params.mask_primer_m1s_z_start_5_Mask_M1S_and_Z_primers.skiprc
//* @style @condition:{method="score",primer_start}, {method="align",maxlen,skiprc} @multicolumn:{method,maxerror,revpr,failed,pf,nrpoc}, {maxlen,skiprc} 

readArray = reads.toString().split(' ')	
revpr_arg = (revpr=="true") ? "--revpr" : ""
failed_arg = (failed=="true") ? "--failed" : "" 

if(method=="score"){
	start = "--start ${primer_start}"
}else{
	start = ""
}

if(method=="align"){
	skiprc_arg = (skiprc=="true") ? "--skiprc" : ""
	maxlen_arg = "--maxlen ${maxlen}"
}else{
	skiprc_arg = ""
	maxlen_arg = ""
}

pf = "--pf ${pf}"

if(mate=="pair"){
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	R1_name =  R1 - '.fastq' - '.fastq'
	R2_name =  R2 - '.fastq' - '.fastq'
	"""
	MaskPrimers.py ${method} -s $R1 -p ${primers} ${start} --mode cut --maxerror ${maxerror} --outname ${R1_name} --log MP_${R1_name}.log  --nproc ${nproc} ${revpr_arg} ${failed_arg} ${skiprc_arg} ${maxlen_arg} ${pf}
	MaskPrimers.py ${method} -s $R2 -p ${primers} ${start} --mode cut --maxerror ${maxerror} --outname ${R2_name} --log MP_${R2_name}.log  --nproc ${nproc} ${revpr_arg} ${failed_arg} ${skiprc_arg} ${maxlen_arg} ${pf}
	"""
}else{
	R1 = readArray[0]
	"""
	echo -e "Assuming inputs for R1"
	MaskPrimers.py ${method} -s $R1 -p ${primers} ${start} --mode cut --maxerror ${maxerror} --outname ${name}_R1 --log MP_R1_${name}.log  --nproc ${nproc} ${revpr_arg} ${failed_arg} ${skiprc_arg} ${maxlen_arg} ${pf}
	"""
}



}


process mask_primer_m1s_z_start_6_Mask_M1S_and_Z_primers {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_primers-pass.fastq$/) "m1s_z_start_6_pass/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /MP_.*$/) "MP_log_start_6/$filename"}
input:
 set val(name),  file(reads) from g39_0_reads_failed2_g1_0
 val mate from g_6_mate_g1_0
 file primers from g0_0_fasta_file0_g1_0

output:
 set val(name), file("*_primers-pass.fastq")  into g1_0_reads0_g_23
 file "MP_*"  into g1_0_logFile1_g1_7
 set val(name), file("*_primers-fail.fastq") optional true  into g1_0_reads_failed2_g13_0

script:
method = params.mask_primer_m1s_z_start_6_Mask_M1S_and_Z_primers.method
maxerror = params.mask_primer_m1s_z_start_6_Mask_M1S_and_Z_primers.maxerror
revpr = params.mask_primer_m1s_z_start_6_Mask_M1S_and_Z_primers.revpr
failed = params.mask_primer_m1s_z_start_6_Mask_M1S_and_Z_primers.failed
pf = params.mask_primer_m1s_z_start_6_Mask_M1S_and_Z_primers.pf
nproc = params.mask_primer_m1s_z_start_6_Mask_M1S_and_Z_primers.nproc
primer_start = params.mask_primer_m1s_z_start_6_Mask_M1S_and_Z_primers.primer_start
maxlen = params.mask_primer_m1s_z_start_6_Mask_M1S_and_Z_primers.maxlen
skiprc = params.mask_primer_m1s_z_start_6_Mask_M1S_and_Z_primers.skiprc
//* @style @condition:{method="score",primer_start}, {method="align",maxlen,skiprc} @multicolumn:{method,maxerror,revpr,failed,pf,nrpoc}, {maxlen,skiprc} 

readArray = reads.toString().split(' ')	
revpr_arg = (revpr=="true") ? "--revpr" : ""
failed_arg = (failed=="true") ? "--failed" : "" 

if(method=="score"){
	start = "--start ${primer_start}"
}else{
	start = ""
}

if(method=="align"){
	skiprc_arg = (skiprc=="true") ? "--skiprc" : ""
	maxlen_arg = "--maxlen ${maxlen}"
}else{
	skiprc_arg = ""
	maxlen_arg = ""
}

pf = "--pf ${pf}"

if(mate=="pair"){
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	R1_name =  R1 - '.fastq' - '.fastq'
	R2_name =  R2 - '.fastq' - '.fastq'
	"""
	MaskPrimers.py ${method} -s $R1 -p ${primers} ${start} --mode cut --maxerror ${maxerror} --outname ${R1_name} --log MP_${R1_name}.log  --nproc ${nproc} ${revpr_arg} ${failed_arg} ${skiprc_arg} ${maxlen_arg} ${pf}
	MaskPrimers.py ${method} -s $R2 -p ${primers} ${start} --mode cut --maxerror ${maxerror} --outname ${R2_name} --log MP_${R2_name}.log  --nproc ${nproc} ${revpr_arg} ${failed_arg} ${skiprc_arg} ${maxlen_arg} ${pf}
	"""
}else{
	R1 = readArray[0]
	"""
	echo -e "Assuming inputs for R1"
	MaskPrimers.py ${method} -s $R1 -p ${primers} ${start} --mode cut --maxerror ${maxerror} --outname ${name}_R1 --log MP_R1_${name}.log  --nproc ${nproc} ${revpr_arg} ${failed_arg} ${skiprc_arg} ${maxlen_arg} ${pf}
	"""
}



}


process mask_primer_m1s_z_align_Mask_M1S_and_Z_primers {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_primers-pass.fastq$/) "m1s_z_align_pass/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /MP_.*$/) "MP_log_align/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_primers-fail.fastq$/) "m1s_z_align_failed/$filename"}
input:
 set val(name),  file(reads) from g1_0_reads_failed2_g13_0
 val mate from g_6_mate_g13_0
 file primers from g0_0_fasta_file0_g13_0

output:
 set val(name), file("*_primers-pass.fastq")  into g13_0_reads0_g_23
 file "MP_*"  into g13_0_logFile1_g13_7
 set val(name), file("*_primers-fail.fastq") optional true  into g13_0_reads_failed22

script:
method = params.mask_primer_m1s_z_align_Mask_M1S_and_Z_primers.method
maxerror = params.mask_primer_m1s_z_align_Mask_M1S_and_Z_primers.maxerror
revpr = params.mask_primer_m1s_z_align_Mask_M1S_and_Z_primers.revpr
failed = params.mask_primer_m1s_z_align_Mask_M1S_and_Z_primers.failed
pf = params.mask_primer_m1s_z_align_Mask_M1S_and_Z_primers.pf
nproc = params.mask_primer_m1s_z_align_Mask_M1S_and_Z_primers.nproc
primer_start = params.mask_primer_m1s_z_align_Mask_M1S_and_Z_primers.primer_start
maxlen = params.mask_primer_m1s_z_align_Mask_M1S_and_Z_primers.maxlen
skiprc = params.mask_primer_m1s_z_align_Mask_M1S_and_Z_primers.skiprc
//* @style @condition:{method="score",primer_start}, {method="align",maxlen,skiprc} @multicolumn:{method,maxerror,revpr,failed,pf,nrpoc}, {maxlen,skiprc} 

readArray = reads.toString().split(' ')	
revpr_arg = (revpr=="true") ? "--revpr" : ""
failed_arg = (failed=="true") ? "--failed" : "" 

if(method=="score"){
	start = "--start ${primer_start}"
}else{
	start = ""
}

if(method=="align"){
	skiprc_arg = (skiprc=="true") ? "--skiprc" : ""
	maxlen_arg = "--maxlen ${maxlen}"
}else{
	skiprc_arg = ""
	maxlen_arg = ""
}

pf = "--pf ${pf}"

if(mate=="pair"){
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	R1_name =  R1 - '.fastq' - '.fastq'
	R2_name =  R2 - '.fastq' - '.fastq'
	"""
	MaskPrimers.py ${method} -s $R1 -p ${primers} ${start} --mode cut --maxerror ${maxerror} --outname ${R1_name} --log MP_${R1_name}.log  --nproc ${nproc} ${revpr_arg} ${failed_arg} ${skiprc_arg} ${maxlen_arg} ${pf}
	MaskPrimers.py ${method} -s $R2 -p ${primers} ${start} --mode cut --maxerror ${maxerror} --outname ${R2_name} --log MP_${R2_name}.log  --nproc ${nproc} ${revpr_arg} ${failed_arg} ${skiprc_arg} ${maxlen_arg} ${pf}
	"""
}else{
	R1 = readArray[0]
	"""
	echo -e "Assuming inputs for R1"
	MaskPrimers.py ${method} -s $R1 -p ${primers} ${start} --mode cut --maxerror ${maxerror} --outname ${name}_R1 --log MP_R1_${name}.log  --nproc ${nproc} ${revpr_arg} ${failed_arg} ${skiprc_arg} ${maxlen_arg} ${pf}
	"""
}



}


process mask_primer_m1s_z_align_parse_log_MP {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tab$/) "m1s_z_log_align_tab/$filename"}
input:
 val mate from g_6_mate_g13_7
 file log_file from g13_0_logFile1_g13_7

output:
 file "*.tab"  into g13_7_logFile00

script:
readArray = log_file.toString()	

"""
ParseLog.py -l ${readArray}  -f ID PRIMER BARCODE ERROR
"""

}


process mask_primer_m1s_z_start_6_parse_log_MP {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tab$/) "m1s_z_log_start_6_tab/$filename"}
input:
 val mate from g_6_mate_g1_7
 file log_file from g1_0_logFile1_g1_7

output:
 file "*.tab"  into g1_7_logFile00

script:
readArray = log_file.toString()	

"""
ParseLog.py -l ${readArray}  -f ID PRIMER BARCODE ERROR
"""

}


process create_m1s_z_reads {

input:
 set val(name), file(read1) from g38_0_reads0_g_23
 set val(name), file(read2) from g39_0_reads0_g_23
 set val(name), file(read3) from g1_0_reads0_g_23
 set val(name), file(read4) from g13_0_reads0_g_23

output:
 set val("M1S"), file("*_M1S.fastq")  into g_23_reads_m1s0_g_28, g_23_reads_m1s0_g35_5
 set val("Z"), file("*_Z.fastq")  into g_23_reads_z1_g_28

script:
reads = read1 + read2 + read3 + read4
read_1_R1 = reads.grep(~/.*R1.*/)[0]
read_1_R2 = reads.grep(~/.*R2.*/)[0]
"""
grep -h -A 3 --no-group-separator  'M1S'  ${read_1_R1} > R1_M1S.fastq
grep -h -A 3 --no-group-separator  'M1S'  ${read_1_R2} > R2_M1S.fastq
grep -h -A 3 --no-group-separator  'Z'  ${read_1_R1} > R1_Z.fastq
grep -h -A 3 --no-group-separator  'Z'  ${read_1_R2} > R2_Z.fastq
"""

}



process Mask_Primer_mask_primers {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_primers-pass.fastq$/) "reads_m1s_umi/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /MP_.*$/) "MP_log_m1s_umi/$filename"}
input:
 set val(name),  file(reads) from g_23_reads_m1s0_g35_5
 file cprimers from g_30_cprimers_g35_5
 file vprimers from g_25_vprimers_g35_5
 val mate from g_6_mate_g35_5

output:
 set val(name), file("*_primers-pass.fastq")  into g35_5_reads00
 file "MP_*"  into g35_5_logFile1_g35_9
 set val(name), file("*_primers-fail.fastq") optional true  into g35_5_reads_failed22

script:
method = params.Mask_Primer_mask_primers.method
maxerror_cprimer = params.Mask_Primer_mask_primers.maxerror_cprimer
maxerror_vprimer = params.Mask_Primer_mask_primers.maxerror_vprimer
revpr = params.Mask_Primer_mask_primers.revpr
failed = params.Mask_Primer_mask_primers.failed
pf = params.Mask_Primer_mask_primers.pf
nproc = params.Mask_Primer_mask_primers.nproc
cprimer_position = params.Mask_Primer_mask_primers.cprimer_position
umi_position = params.Mask_Primer_mask_primers.umi_position
barcode_both = params.Mask_Primer_mask_primers.barcode_both
cprimer_start = params.Mask_Primer_mask_primers.cprimer_start
vprimer_start = params.Mask_Primer_mask_primers.vprimer_start
umi_length = params.Mask_Primer_mask_primers.umi_length
maxlen = params.Mask_Primer_mask_primers.maxlen
skiprc = params.Mask_Primer_mask_primers.skiprc
//* @style @condition:{method="score",cprimer_start,vprimer_start,umi_length}, {method="align",maxlen,skiprc} @multicolumn:{method,maxerror_cprimer,maxerror_vprimer,revpr,failed,pf,nrpoc},{cprimer_position,umi_position,barcode_both}, {cprimer_start,vprimer_start,umi_length}, {maxlen,skiprc} 

readArray = reads.toString().split(' ')	
revpr_arg = (revpr=="true") ? "--revpr" : ""
failed_arg = (failed=="true") ? "--failed" : "" 
pf = "--pf ${pf}"

R1_barcode = ""
R2_barcode = ""

R1_maxerror = 0.2
R2_maxerror = 0.2
if(cprimer_position=='R1'){
	R1_primers = "${cprimers}"
	R2_primers = "${vprimers}"
	
	R1_maxerror = maxerror_cprimer
	R2_maxerror = maxerror_vprimer
	
	if(umi_position == 'R1'){
		R1_start = umi_length + cprimer_start
		R1_barcode = "--barcode"
		R2_start = vprimer_start
	}else{
		R1_start = cprimer_start
		R2_start = vprimer_start + umi_length
		R2_barcode = "--barcode"
	}
}else{
	R1_primers = "${vprimers}"
	R2_primers = "${cprimers}"
	
	R2_maxerror = maxerror_cprimer
	R1_maxerror = maxerror_vprimer
	
	if(umi_position == 'R1'){
		R1_start = umi_length + vprimer_start
		R1_barcode = "--barcode"
		R2_start = cprimer_start
	}else{
		R1_start = vprimer_start
		R2_start = cprimer_start + umi_length
		R2_barcode = "--barcode"
	}
}

if(barcode_both=="true"){
	R1_barcode = "--barcode"
	R2_barcode = "--barcode"
	if(umi_position == 'R1'){
		R2_start = R2_start + umi_length
	}else{
		R1_start = R1_start + umi_length
	}
}

R1_start = (method=="align") ? "" : "--start ${R1_start}"
R2_start = (method=="align") ? "" : "--start ${R2_start}"

if(method=="align"){
	skiprc_arg = (skiprc=="true") ? "--skiprc" : ""
	maxlen_arg = "--maxlen ${maxlen}"
}else{
	skiprc_arg = ""
	maxlen_arg = ""
}

if(mate=="pair"){
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	
	"""
	MaskPrimers.py ${method} -s $R1 -p ${R1_primers} ${R1_start} ${R1_barcode} --mode cut --maxerror ${R1_maxerror} --outname ${name}_R1 --log MP_R1_${name}.log  --nproc ${nproc} ${revpr_arg} ${failed_arg} ${skiprc_arg} ${maxlen_arg} ${pf}
	MaskPrimers.py ${method} -s $R2 -p ${R2_primers} ${R2_start} ${R2_barcode} --mode cut --maxerror ${R2_maxerror} --outname ${name}_R2 --log MP_R2_${name}.log  --nproc ${nproc} ${revpr_arg} ${failed_arg} ${maxlen_arg} ${skiprc_arg} ${pf}
	"""
}else{
	R1 = readArray[0]
	R1_primers = R1_primers - '.fastq' - '.fastq'
	R1_primers = '${R1_primers}.fasta'
	"""
	echo -e "Assuming inputs for R1"
	MaskPrimers.py ${method} -s $R1 -p ${R1_primers} ${R1_start} ${R1_barcode} --mode cut --maxerror ${R1_maxerror} --outname ${name}_R1 --log MP_R1_${name}.log  --nproc ${nproc} ${revpr_arg} ${failed_arg} ${maxlen_arg} ${skiprc_arg} ${pf}
	"""
}



}


process Mask_Primer_parse_log_MP {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tab$/) "MP_umi_log_tab/$filename"}
input:
 val mate from g_6_mate_g35_9
 file log_file from g35_5_logFile1_g35_9

output:
 file "*.tab"  into g35_9_logFile00

script:
readArray = log_file.toString()	

"""
ParseLog.py -l ${readArray}  -f ID PRIMER BARCODE ERROR
"""

}


process PairAwk {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /reads\/.*q$/) "reads_m1s_z/$filename"}
input:
 set val(name_1), file(reads_1) from g_23_reads_m1s0_g_28
 set val(name_2), file(reads_2) from g_23_reads_z1_g_28
 val mate from g_6_mate_g_28
 file pairawk from g_29_PAIRAWK_script_g_28

output:
 set val("M1S_Z"), file("reads/*q")  into g_28_reads0_g_46

script:

if(mate=="pair"){
	read_1_R1 = reads_1.grep(~/.*R1.*/)[0]
	read_1_R2 = reads_1.grep(~/.*R2.*/)[0]
	read_2_R1 = reads_2.grep(~/.*R1.*/)[0]
	read_2_R2 = reads_2.grep(~/.*R2.*/)[0]
	"""
	./${pairawk} ${read_1_R1} ${read_2_R2}
	./${pairawk} ${read_1_R2} ${read_2_R1}
	
	mkdir -p reads
	
	ls ./*_pair-pass* | grep "M1S"
	ls ./*_pair-pass* | grep "M1S" | xargs cat > reads/R1_pair-pass.fastq
	ls ./*_pair-pass* | grep "Z"
	ls ./*_pair-pass* | grep "Z" | xargs cat > reads/R2_pair-pass.fastq  			
	"""
}else{
	"""
	echo -e 'PairAwk works only on pair-end reads.'
	"""
}

}


process split_seq_samples {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /samples\/.*$/) "samples/$filename"}
input:
 set val(name),file(reads) from g_28_reads0_g_46

output:
 set val(name), file("samples/*")  into g_46_reads00

script:
field = params.split_seq_samples.field

reads = reads.join(" ")

"""
#!/bin/bash

mkdir -p samples
SplitSeq.py group -s ${reads} -f ${field} --outdir samples

for file in samples/*
do 
	filename=\${file##*/}
	filename=\${filename%.*}
	part1=\${filename%_"${field}"*}
	part2=\${filename#*_"${field}"-}
	
	echo \$part2"_"\$part1".fastq"
	
	mv \$file samples/\$part2"_"\$part1".fastq"
done
"""

}


process mask_primer_m1s_z_start_5_parse_log_MP {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tab$/) "m1s_z_log_start_5_tab/$filename"}
input:
 val mate from g_6_mate_g39_7
 file log_file from g39_0_logFile1_g39_7

output:
 file "*.tab"  into g39_7_logFile00

script:
readArray = log_file.toString()	

"""
ParseLog.py -l ${readArray}  -f ID PRIMER BARCODE ERROR
"""

}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
