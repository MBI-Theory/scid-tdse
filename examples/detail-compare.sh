#!/bin/bash
sum="$1"
ref="${2:-${sum}_ref}"
#
gnuplot <<-eoi
        set terminal pdf
	set output "detail-compare.pdf"
	set grid
	do for [c=3:32] {
	  if (c== 3) { set title "Vector-potential magnitude" }
	  if (c== 4) { set title "Vector-potential, lab theta" }
	  if (c== 5) { set title "Vector-potential, lab phi" }
	  if (c== 6) { set title "Re[<Psi_L|Psi_R>]" }
	  if (c== 7) { set title "Im[<Psi_L|Psi_R>]" }
	  if (c== 8) { set title "Re[<Psi_L|H_at+H_L+V_cap|Psi_R>], Hartree" }
	  if (c== 9) { set title "Im[<Psi_L|H_at+H_L+V_cap|Psi_R>], Hartree" }
	  if (c==10) { set title "Re[<Psi_L|H_at+H_L|Psi_R>], Hartree" }
	  if (c==11) { set title "Im[<Psi_L|H_at+H_L|Psi_R>], Hartree" }
	  if (c==12) { set title "Re[<Psi_L|e . x|Psi_R>], e-Bohr" }
	  if (c==13) { set title "Im[<Psi_L|e . x|Psi_R>], e-Bohr" }
	  if (c==14) { set title "Re[<Psi_L|e . y|Psi_R>], e-Bohr" }
	  if (c==15) { set title "Im[<Psi_L|e . y|Psi_R>], e-Bohr" }
	  if (c==16) { set title "Re[<Psi_L|e . z|Psi_R>], e-Bohr" }
	  if (c==17) { set title "Im[<Psi_L|e . z|Psi_R>], e-Bohr" }
	  if (c==18) { set title "Re[(d^2/d t^2) <Psi_L|e . x|Psi_R>], e-Bohr/jiffy^2" }
	  if (c==19) { set title "Im[(d^2/d t^2) <Psi_L|e . x|Psi_R>], e-Bohr/jiffy^2" }
	  if (c==20) { set title "Re[(d^2/d t^2) <Psi_L|e . y|Psi_R>], e-Bohr/jiffy^2" }
	  if (c==21) { set title "Im[(d^2/d t^2) <Psi_L|e . y|Psi_R>], e-Bohr/jiffy^2" }
	  if (c==22) { set title "Re[(d^2/d t^2) <Psi_L|e . z|Psi_R>], e-Bohr/jiffy^2" }
	  if (c==23) { set title "Im[(d^2/d t^2) <Psi_L|e . z|Psi_R>], e-Bohr/jiffy^2" }
	  if (c==24) { set title "Re[(d/d t) <Psi_L|e . x|Psi_R>], e-Bohr/jiffy" }
	  if (c==25) { set title "Im[(d/d t) <Psi_L|e . x|Psi_R>], e-Bohr/jiffy" }
	  if (c==26) { set title "Re[(d/d t) <Psi_L|e . y|Psi_R>], e-Bohr/jiffy" }
	  if (c==27) { set title "Im[(d/d t) <Psi_L|e . y|Psi_R>], e-Bohr/jiffy" }
	  if (c==28) { set title "Re[(d/d t) <Psi_L|e . z|Psi_R>], e-Bohr/jiffy" }
	  if (c==29) { set title "Im[(d/d t) <Psi_L|e . z|Psi_R>], e-Bohr/jiffy" }
	  if (c==30) { set title "Electric field, lab X, atomic units" }
	  if (c==31) { set title "Electric field, lab Y, atomic units" }
	  if (c==32) { set title "Electric field, lab Z, atomic units" }
	  plot "$ref" u 2:c w l lw 1.5, "$sum" u 2:c w l lw 0.5
	  }
eoi
