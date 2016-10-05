BEGIN {
    PROCINFO["sorted_in"]="@ind_str_asc" # traversal order for for(i in a)
    if(cols) {                           # set flag -v cols="bxbam_columns.txt"
        while ((getline line< cols)>0) { # parse line by line
            gsub(/: [^ ]+/,"",line)      # remove values from "key: value" 
            split(line,a)                # split to temporary array  ## consider parsing key:value:value
            for(i in a)                  # collect keys to column array
                col[a[i]]
        }
        for(i in col)                    # output columns
            printf "%6s%s", i, OFS
        print ""
    }
}
NR==1 {
    next
}
{
    gsub(/: /,"=")                       # replaces key-value separator ": " with "="
    split($0,b,FS)                       # splits record from separator FS
    for(i in b) {
        split(b[i],c,"=")                # splits key=value to c[1]=key, c[2]=value
        b[c[1]]=c[2]                     # b[key]=value
    }
    for(i in col)                        # go throurgh headers in col[] and printf from b[]
        printf "%6s%s", (i in b?b[i]:"NaN"), OFS; print ""
}

# awk -v cols="bxbam_columns.txt" -f bxbam.awk bxbam_columns.txt

# NR =  total number of records processed
# FNR =  total number of records for each input file, typically the line number
