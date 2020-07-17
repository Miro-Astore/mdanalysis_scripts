 set ref_sel [atomselect top "name CA" frame 0 ]
puts hi
 set num_frames [molinfo top get numframes]

 for {set i 0} {$i < $num_frames} {incr i} {
		 animate goto $i
		 set m_sel [atomselect top "name CA"]
		 set a_sel [atomselect top "all"]
		 set M [measure fit $m_sel $ref_sel]
		 $a_sel move $M
}
