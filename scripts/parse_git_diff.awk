#!/usr/bin/awk -f

/^diff --git/ {
    file = $4
    sub(/^b\//, "", file)
}

/^@@.*\+([0-9]+)(,([0-9]+))?/ {
    match($0, /\+([0-9]+)(,([0-9]+))?/, arr)
    a = arr[1]
    l = (arr[3] != "") ? arr[3] : 1
    for (i = 0; i < l; i++) {
        print "file=" file ",line=" a + i ","
    }
}
