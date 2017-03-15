inName=$1

n=1
for file in `ls ${inName}*.png`; do
    echo $file
    cp $file links/imLink_`printf "%02d" ${n}`.png; n=$((n+1))
done

cd links
ffmpeg -r 5 -f image2 -i imLink_%02d.png ../${inName}.mp4
rm -f imLink*
cd ../

