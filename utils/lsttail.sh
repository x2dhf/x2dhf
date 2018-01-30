#! /bin/bash

################################################################################
#                                                                              #
#  Copyright (C) 2010 Jacek Kobus <jkob@fizyka.umk.pl>                         #
#                                                                              #
#  This software may be used and distributed according to the terms            #
#  of the GNU General Public License, see README and COPYING.                  #
#                                                                              #
################################################################################

# By default this script prints the last 100 lines of the newest *.lst
# file found in a current directory

while getopts ":f:hn:" Option
do
  case $Option in
    f ) f=$OPTARG;;
    h ) echo "Usage: lsttail.sh [-f <file>] [-h] [-n <number of lines>]";exit ;;
    n ) n=$OPTARG;;
  esac
done

num=${n:-100}
file=${f:-`ls *.lst`}

newest_lst=`ls -t $file | head -1`

tail -$num $newest_lst

# 




