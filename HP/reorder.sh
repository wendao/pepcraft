L=$1
awk '
NR==FNR {order[$1] = NR; next}
{
    if ($1 in order) {
        print order[$1], $0
    } else {
        print length(order)+1, $0
    }
}
' rand_seq${L}.txt rand_seq${L}.out.txt | sort -n -k1 | cut -d' ' -f2- > L${L}_sorted

#awk '$3==0{print $2}' rand_seq${L}.txt | sort | uniq -c | sort -nk2 > E_L${L}_undes.hist
#awk '$3==1{print $2}' rand_seq${L}.txt | sort | uniq -c | sort -nk2 > E_L${L}_des.hist

