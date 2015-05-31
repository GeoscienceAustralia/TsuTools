## Various comparisons between SFFM with normal phase-generating-field, and the
## corresponding SFFM which uses the stable distribution

normal_dist_runs = Sys.glob('OUTPUTS_2/Rimages/*_S_G*')
stable_dist_runs = Sys.glob('OUTPUTS_2/Rimages/*_S_S*')

model_labels=c('A_', 'C_', 'AF', 'CF')

for(i in 1:length(model_labels)){

    normal_run = normal_dist_runs[grep(model_labels[i], normal_dist_runs)]
    stable_run = stable_dist_runs[grep(model_labels[i], stable_dist_runs)]

    print('Comparing:')
    print(normal_run)
    print(stable_run)

    eN = new.env()
    eS = new.env()

    load(normal_run, envir = eN)
    load(stable_run, envir = eS)

    figure_title=paste0('alpha_output_results/stable_vs_normal_median80plus_', model_labels[i], '.pdf')

    pdf(figure_title, width=10, height=6)

    # Compare median 80+ slip and q95 80+ slip for synthetic runs
    for(j in 1:66){
        par(mfrow=c(1,2))
        qqplot(eN$rod$myrast_simstats[[j]]$median_80_slip,
               eS$rod$myrast_simstats[[j]]$median_80_slip,
               xlab='Normal', ylab = 'Stable')
        abline(0,1,col='red')
        qqplot(eN$rod$myrast_simstats[[j]]$quant_95_80_slip,
               eS$rod$myrast_simstats[[j]]$quant_95_80_slip,
               xlab='Normal', ylab = 'Stable')
        abline(0,1,col='red')
    } 

    dev.off()

    figure_title=paste0('alpha_output_results/Corner_wavenumbers_stable_vs_normal_', model_labels[i], '.pdf')
    pdf(figure_title, width=10, height=10)
    par(mfrow=c(2,2))
    plot(eN$rod$p1$kcxN, eS$rod$p1$kcxN, main='Numerical kcx')
    abline(0,1,col='red')
    plot(eN$rod$p1$kcyN, eS$rod$p1$kcyN, main='Numerical kcy')
    abline(0,1,col='red')

    plot(eN$rod$p1$kcx, eS$rod$p1$kcx, main='Physical kcx')
    abline(0,1,col='red')
    plot(eN$rod$p1$kcy, eS$rod$p1$kcy, main='Physical kcy')
    abline(0,1,col='red')
    dev.off() 
}


