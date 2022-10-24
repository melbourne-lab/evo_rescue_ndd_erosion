### Combine simulations from September 22.

all.data = 
  rbind(
    read.csv('simulations/outputs/alldata_n100_a000_hivar.csv'),
    read.csv('simulations/outputs/alldata_n100_a000_lowvar.csv'),
    read.csv('simulations/outputs/alldata_n20_a000_hivar.csv'),
    read.csv('simulations/outputs/alldata_n20_a000_lowvar.csv'),
    read.csv('simulations/outputs/alldata_n100_a035_hivar.csv'),
    read.csv('simulations/outputs/alldata_n100_a035_lowvar.csv'),
    read.csv('simulations/outputs/alldata_n20_a035_hivar.csv'),
    read.csv('simulations/outputs/alldata_n20_a035_lowvar.csv')
  )

write.csv(all.data, row.names = FALSE,
          'simulations/outputs/alldata_combined.csv')
