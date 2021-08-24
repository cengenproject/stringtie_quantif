BEGIN{
	FS=OFS="\t"
}

{
	if($7~"\\."){
		# when strand is not defined ("."), duplicate line and make a + and - line.
		gsub("\\.","+",$7)
		$9=gensub("(MSTRG[0-9.]+)", "\\1.f", "g", $9)
		print $0
		gsub("\+","-",$7)
		$9=gensub("(MSTRG[0-9.]+)\\.f", "\\1.r", "g",$9)
		print $0
	}else {
		if($2~"rtracklayer"){
			
			$9=gensub("transcript_id \"(mRNA[0-9]+)\"", "gene_id \"rtck_\\1\"; transcript_id \"\\1\"", "g", $9)
		}
		print $0
	}
}


