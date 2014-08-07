 set term post eps color solid
 set term post enhanced
 set term post "Times-Roman" 22
#set  term post "Times-Roman" 20 fontfile "/home/youjin/fonts/cmsy10.pfa"
#set output "qn.eps"
#set term post"Times-Roman" 12
#set output "Sigma_compare.eps"
set output "Re_Sigma.eps"
set xlabel "tau" 
set ylabel "Re(Sigma)"
#set title ""
#set zlabel " L"
#set logscale 
#set logscale x
#set key graph 0.85, graph 0.85
#set xtics (0.0084275,0.0084335,0.0084395)
#set xrange [0:127]
#set size square
#set label "F" at 0.2, 0.2 font "Times-Roman,15"
#set nokey
set pointsize 0.5
set key right top
#set key right bottom 
#set key left bottom box 0 samplen 0.2
#set key left bottom box 0 samplen 0.645,0.1
plot "0.50_Sigma_0.dat" us 1:(-1.0*$2) notit wi poi lt 1, \
		 "0.50_Sigma_0.dat" us 1:(-1.0*$2) tit "0 mc" wi l lt 1, \
     "0.50_Sigma_1_full.dat" us 1:(-1.0*$2) notit wi poi lt 3, \
		 "0.50_Sigma_1_full.dat" us 1:(-1.0*$2) tit "1 mc" wi l lt 3, \
     "0.50_Sigma_2_full.dat" us 1:(-1.0*$2) notit wi poi lt 5, \
		 "0.50_Sigma_2_full.dat" us 1:(-1.0*$2) tit "2 mc" wi l lt 5, \
     "0.50_Sigma_3_full.dat" us 1:(-1.0*$2) notit wi poi lt 7, \
		 "0.50_Sigma_3_full.dat" us 1:(-1.0*$2) tit "3 mc" wi l lt 7, \
		 "Sigma_0.dat" us 1:2 notit wi poi lt 2, "Sigma_0.dat" us 1:2 tit "0 Nikolay" wi l lt 2, \
		 "Sigma_1_full.dat" us 1:2 notit wi poi lt 4, "Sigma_1_full.dat" us 1:2 tit "1 Nikolay" wi l lt 4, \
		 "Sigma_2_full.dat" us 1:2 notit wi poi lt 6, "Sigma_2_full.dat" us 1:2 tit "2 Nikolay" wi l lt 6, \
		 "Sigma_3_full.dat" us 1:2 notit wi poi lt 8, "Sigma_3_full.dat" us 1:2 tit "3 Nikolay" wi l lt 8
#smooth csplines notitle  w linespoNikolays
pause -1
#plot "0.50_Gam_MC.dat" us 1:3:($3*$4) notit wi err lt 1,"0.50_Gam_MC.dat" us 1:3 tit "1 mc" wi l lt 1, \