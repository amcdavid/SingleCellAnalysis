setClass('PairwisePearson', contains='matrix', representation=representation(nsamp='numeric', npos='numeric'), validity=function(object){
    ncol(object@.Data)==nrow(object@.Data) && ncol(object@.Data) == length(object@npos)
})

setClass('PairwiseOddsratio', contains='matrix', representation=representation(genes='character'))

setClass('PairwiseOddsWithHook', contains='PairwiseOddsratio', representation=representation(hook='matrix'))
