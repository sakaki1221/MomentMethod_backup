
      set format x "%.2e"

      set xlabel "a1"
      set ylabel "erg"
      plot "./k.txt" using 1:2 title "k" w lp

      replot "./k.txt" using 1:3 title "term 1,    2*de2dr2(a1)" w lp

      replot "./k.txt" using 1:4 title "term 2,    4*dedr(a1)/a1" w lp

      replot "./k.txt" using 1:5 title "term 3,    de2dr2(a2)" w lp

      replot "./k.txt" using 1:6 title "term 4,    2*dedr(a2)/a2" w lp

      