echo "print '.' * 100" | python | RNAblueprint -s 1 -n 10  | while read sequence; do RNAfold | grep '(' | cut -f1 -d" "; done > structures.txt
for i in $(< structures.txt); do echo $i | RNAblueprint -s 1 -n 100; done > results.txt
