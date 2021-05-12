#测试电脑相对性能
library(benchmarkme)
res_io = benchmark_io(runs = 1, size = 5)
upload_results(res_io)
plot(res_io)

#测试CPU性能
res = benchmark_std(runs = 1)
upload_results(res)
plot(res)

