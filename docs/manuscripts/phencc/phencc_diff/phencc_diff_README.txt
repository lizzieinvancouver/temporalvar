Started 30 April 2021

I pulled figures as submitted previously if I thought they might change and put in phencc_diff/figs_maintext. So far they have not changed so I just added ..// to the paths. 

I ended up putting the .tex files all in one folder, because otherwise the paths were annoying. Anyway, below is the code and steps. 

https://texblog.org/2018/08/14/track-changes-with-latexdiff/

(1) Put the new files into phencc_diff and add ..// as needed to file paths. Check that both the new file and the phencc_old.tex compile on their own okay. 

(2) Do this in Terminal: 

cd projects/temporalvar/docs/manuscripts/phencc/phencc_diff

latexdiff phencc_old.tex phencc.tex > phencc_diff.tex

(3) You can then just open phencc_diff.tex and compile, bib text it, compile it again in Aquamacs or some other latex programs. 

* This all worked well on the first round but was upset on the second -- it somehow really did see an issue on line 128 that it stopped it compiling. I used texshop (better line numbering) to confirm and ended up pasting in the old code for that tiny chunk from git.* 