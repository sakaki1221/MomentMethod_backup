
      set format x "%.2e"

      set xlabel "a1"
      set ylabel "erg"
      plot "./gamma.txt" using 1:2 title "gamma" w lp

      replot "./gamma.txt" using 1:3 title "xxxx" w lp

      replot "./gamma.txt" using 1:4 title "xxyy" w lp

      