 set term post eps color solid
 set term post enhanced
 set term post "Times-Roman" 22
#set  term post "Times-Roman" 20 fontfile "/home/youjin/fonts/cmsy10.pfa"
#set output "qn.eps"
#set term post"Times-Roman" 12
#set output "Sigma_compare.eps"
set output "Re_Gamma.eps"
set xlabel "tau1=tau2" 
set ylabel "Re(Gamma)"
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
plot "0.90_1_coll_GamMC.dat" us 1:2 notit wi poi lt 1, \
		 "0.90_1_coll_GamMC.dat" us 1:2 tit "1 mc" wi l lt 1, \
		 "0.90_Gam1.dat" us 1:2 notit wi poi lt 3, "0.90_Gam1.dat" us 1:2 tit "1 int" wi l lt 3
#smooth csplines notitle  w linespoints
pause -1
#plot "0.50_1_coll_GamMC.dat" us 1:2 notit wi poi lt 1, \
		 #"0.50_1_coll_GamMC.dat" us 1:2 tit "1 mc" wi l lt 1, \
		 #"0.50_Gam1.dat" us 1:2 notit wi poi lt 3, \
		 #"0.50_Gam1.dat" us 1:2 tit "1 int" wi l lt 3
