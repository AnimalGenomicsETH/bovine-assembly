#pigz -dc $1 |
cat $1 |
     awk 'NR%2==0{c++; l+=length($0)}
          END{
                print "Number of reads: "c; 
                print "Number of bases in reads: "l
              }'
