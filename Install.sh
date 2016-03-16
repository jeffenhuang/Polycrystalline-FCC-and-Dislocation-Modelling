mode=$1

action() {
  if (test $mode=0) then
    rm -f ../$1
  elif (! cmp -s $1 ../$1) then
    if (test -z "$2" || test -e ../$2) then
    	cp $1 ..
    	if (test $mode = 2) then
    	  echo "  updating src/$1"
    	fi
    fi
  elif (test -n "$2") then
    if (test ! -e ../$2) then
      rm -f ../$1
    fi
  fi
}

for file in *.cpp *.h; do
  action $file
done
