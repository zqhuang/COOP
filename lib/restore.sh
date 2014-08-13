#! /bin/bash
for i in `ls camb/*__.bak`
do
    cp ${i} ${i/__.bak/}
    echo "${i/__.bak/} restored"
done
for i in `ls source/*__.bak`
do
    cp ${i} ${i/__.bak/}
    echo "${i/__.bak/} restored"
done
rm -f camb/*__.bak
rm -f source/*__.bak


