pigz -dc $1 |
     awk 'NR%4==2{c++; l+=length($0)}
          END{
                print "Number of reads: "c; 
                print "Number of bases in reads: "l
              }'
