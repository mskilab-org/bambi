  {
      for (i=1; i<=NF; i++) {
          match($i, /^([^:]+):(.*)/, m)
          a[NR,m[1]] = m[2]
          cols[m[1]]
     }
  }
  END{n=asorti(cols,icols);
      for(j=1;j<=n;j++) printf "%s", icols[j] OFS;
      print "";
      for(i=1;i<=NR;i++) # print
          {for(j=1;j<=n;j++)
              {v=a[i,icols[j]];
              printf "%s", (v?v:"NaN") OFS}
           print ""}
  }
