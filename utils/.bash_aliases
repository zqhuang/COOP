alias em='emacs --no-splash'
alias open="evince"
alias hm='cd ~'
alias wk='cd ~/work'
alias rm='rm -I'
alias got='cd ~/work/talks/beamer/'
alias coop='cd ~/work/scilibs/COOP'
alias gw='ssh gw.cita.utoronto.ca'
alias svngw='ssh -L 5902:kodiak.cita.utoronto.ca:22 gw.cita.utoronto.ca'
alias gpull='git pull'
alias gci='git commit'
alias gco='git checkout'
alias gpush='git push'
alias gstat='git status'
alias cori='ssh -l zqhuang cori.nersc.gov'
alias edison='ssh -l zqhuang edison.nersc.gov'
alias scinet='ssh -l zqhuang login.scinet.utoronto.ca'
alias magique='ssh -l zhuang01 magique3.iap.fr'
alias fasy='/home/zqhuang/work/scilibs/COOP/utils/fasy.sh'
alias findp='/home/zqhuang/work/papers/search.py'
alias editp='/home/zqhuang/work/papers/maintain.py'
export CITA='gw.cita.utoronto.ca'
export CLICOLOR=1
export LS_COLORS='di=0;35:ln=36:or=31:mi=31:ex=32'
export PATH=$PATH:/home/zqhuang/bin/
export CORI='zqhuang@cori.nersc.gov'
export EDISON='zqhuang@edison.nersc.gov'
export MAGIQUE='zhuang01@magique3.iap.fr'
export SCINET='zqhuang@login.scinet.utoronto.ca'
export PS1='zq$ '

export MAS_ROOT='/home/zqhuang/work/scilibs/mce_script'

#alias svnco='svn co svn+ssh://myhost/home/zqhuang/work/svnrep/'


# cfitsio built separately from healpix?
export EXTERNAL_CFITSIO=yes

# install prefix (lib and include will be appended as necessary)
export CFITSIO_EXT_PREFIX=/home/zqhuang/work/scilibs/cfitsio/

# healpix install directory
export HEALPIX=/home/zqhuang/work/scilibs/Healpix/