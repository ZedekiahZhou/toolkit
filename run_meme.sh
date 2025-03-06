#/bin/bash
# micromamba activate meme
for file in *.fa
do
  if [ ! -d "${file%.fa}" ]; then
    meme ${file} -o ${file%.fa} -dna -w 5 -p 10
  fi
done
