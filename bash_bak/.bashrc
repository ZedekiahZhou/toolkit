export HOME=/lustre2/chengqiyi_pkuhpc/zhouzhe
#export R_LIBS='~/software/rpkgs'
alias condabio='conda activate zzhou_bio'
export PATH='/lustre2/chengqiyi_pkuhpc/zhouzhe/miniconda3/bin/':$PATH
export PATH='~/toolkit/':$PATH
source /lustre2/chengqiyi_pkuhpc/zhouzhe/toolkit/custom_commands.sh

LANG=en_US.UTF-8

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/lustre2/chengqiyi_pkuhpc/zhouzhe/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/lustre2/chengqiyi_pkuhpc/zhouzhe/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/lustre2/chengqiyi_pkuhpc/zhouzhe/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/lustre2/chengqiyi_pkuhpc/zhouzhe/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup

if [ -f "/lustre2/chengqiyi_pkuhpc/zhouzhe/miniconda3/etc/profile.d/mamba.sh" ]; then
    . "/lustre2/chengqiyi_pkuhpc/zhouzhe/miniconda3/etc/profile.d/mamba.sh"
fi
# <<< conda initialize <<<

condabio
