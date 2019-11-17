#!/bin/sh

# exe_入力ﾌｧｲﾙ_流域面積_出力ﾌｧｲﾙ

# PROSY="2000 2008 2009 2010 2011 2012"
# PROSY="2000 2012"
PROSY="2100"

INDAT="InputData.txt"

#=== 期間を指定し,年別にﾃﾞｰﾀを抽出
  for YYYY in ${PROSY}
  do
    ST=${YYYY}040101
    ET=${YYYY}053124
    FIL=./Input${YYYY}.txt
    echo "FIL-->"${FIL}
    awk '{
      if($(1)>=stim && $(1)<=etim){ print $0; }
    }' stim=${ST} etim=${ET} ${FIL} > ./tmp/dt_${YYYY}.tmp
  done


# ﾌｧｲﾙを結合
  if [ -e ${INDAT} ] ; then
    rm ${INDAT}
    echo ' Delete >>  '${INDAT}
  fi

  for YYYY in ${PROSY}
  do
    LL=`cat ./tmp/dt_${YYYY}.tmp | wc -l`
    # LL=`wc -l ./tmp/dt_${YYYY}.tmp | cut -d' ' -f1`
    echo ${LL} >> ${INDAT}
    cat ./tmp/dt_${YYYY}.tmp >> ${INDAT}
  done
#  echo "-999  -999  -999" >> ${INDAT} # エラーになるので修正
  echo "-999" >> ${INDAT}

#  rm ./tmp/*


#  ./FitTank.exe InputData.txt 134.0 sol_All.csv
  ./FitTank InputData.txt 134.0 sol_All.csv
