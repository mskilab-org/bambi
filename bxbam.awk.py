$ cat > columns.txt  # use three separate columns files for BX, MD&BX, all? Best option---how to use parameter?
BEGIN {
    PROCINFO["sorted_in"]="@ind_str_asc" # traversal order for for(i in a)
    if(cols) {                           # set flag -v cols="bxbam_columns.txt"
        while ((getline line< cols)>0) { # read it in line by line
            gsub(/: [^ ]+/,"",line)      # remove values from "key: value"
            split(line,a)                # split to temp array
            for(i in a)                  # collect keys to column array
                col[a[i]]
        }
        for(i in col)                    # output columns
            printf "%6s%s", i, OFS
        print ""
}
}
NR==1 && cols=="" {                      # if the header cols are in the beginning of data file
    # if not, -v cols="column_file.txt"
    split($0, a, " +")                   # split header record by spaces
    for(i in a) {
        col[a[i]]                        # set them to array col
        printf "%6s%s", a[i], OFS        # output the header
    }
    print ""
}
NR==1 {
    next
}
{
    gsub(/: /,"=")                       # replace key-value separator ": " with "="
    split($0,b,FS)                       # split record from separator FS
    for(i in b) {
        split(b[i],c,"=")                # split key=value to c[1]=key, c[2]=value
        b[c[1]]=c[2]                     # b[key]=value
    }
    for(i in col)                        # go thru headers in col[] and printf from b[]
        printf "%6s%s", (i in b?b[i]:"NaN"), OFS; print ""
}

