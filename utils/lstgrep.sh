#! /bin/bash

################################################################################
#                                                                              #
#  Copyright (C) 2010 Jacek Kobus <jkob@fizyka.umk.pl>                         #
#                                                                              #
#  This software may be used and distributed according to the terms            #
#  of the GNU General Public License, see README and COPYING.                  #
#                                                                              #
################################################################################

# By default this script prints the last 100 lines matching 'total' of the newest *.lst
# file found in a current directory

while getopts ":f:hn:s:" Option
do
  case $Option in
    f ) f=$OPTARG;;
    h ) echo "Usage: lstgrep.sh [-f <file>] [-h] [-n <number of lines>] [-s keyword]";exit ;;
    n ) n=$OPTARG;;
    s ) s="$OPTARG";;
  esac
done

file=${f:-`ls *.lst`}
keyword=${s:-total}
num=${n:-100}

newest_lst=`ls -t $file | head -1`

grep "$keyword" $newest_lst | tail -$num
# 




