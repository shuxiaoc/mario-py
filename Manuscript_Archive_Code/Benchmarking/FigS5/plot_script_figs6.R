## bokai zhu
##

# this script is related to the script making for the figureS6
# figures related to the benchmark of MARIO with time/memory/sparsity

library(ggplot2)
############################ time performance ##################################

### full pipeline time performance

full_time = read.csv("/home/bkzhu/SNE-multi/figure_rcode/sup_simulation/code_0713/data/timeit_full_fig2/timeit_full.csv")
full_time$total = full_time$n1 + full_time$n2
full_time$time = full_time$time / 60
p = ggplot(data=full_time, aes(x=total, y=time)) +
  geom_line()+
  geom_point()+ theme_classic() + ggtitle("MARIO full pipeline Time")
ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/sup_simulation/code_0713/data/figs6/fig2full_time.svg", plot=p, width=6, height=3)

full_time = read.csv("/home/bkzhu/SNE-multi/figure_rcode/sup_simulation/code_0713/data/timeit_full_fig3/timeit_full.csv")
full_time$total = full_time$n1 + full_time$n2
full_time$time = full_time$time / 60
p = ggplot(data=full_time, aes(x=total, y=time)) +
  geom_line()+
  geom_point()+ theme_classic() + ggtitle("MARIO full pipeline Time")
ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/sup_simulation/code_0713/data/figs6/fig3full_time.svg", plot=p, width=6, height=3)

full_time = read.csv("/home/bkzhu/SNE-multi/figure_rcode/sup_simulation/code_0713/data/timeit_full_fig4/timeit_full.csv")
full_time$total = full_time$n1 + full_time$n2
full_time$time = full_time$time / 60
p = ggplot(data=full_time, aes(x=total, y=time)) +
  geom_line()+
  geom_point()+ theme_classic() + ggtitle("MARIO full pipeline Time")
ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/sup_simulation/code_0713/data/figs6/fig4full_time.svg", plot=p, width=6, height=3)


##### only matching steps
## first part the plot for time, fig2
mario_time = read.csv("/home/bkzhu/SNE-multi/figure_rcode/sup_simulation/code_0713/data/timeit_fig2/timeit.csv")
mario_time$totalcell = mario_time$n1 + mario_time$n2
mario_time$ovlp_all = mario_time$time_ovlp + mario_time$time_all
mario_time$sparse_level = rep(c("lv1","lv1","lv2","lv2","lv3","lv3","lv4","lv4","lv5","lv5"),5)
mario_time = subset(mario_time,mario_time$mode == "sparse") # dont need dense time since full sparse is dense
mario_time = subset(mario_time,mario_time$sparse_level %in% c("lv1","lv2","lv5") )# dont need dense time since full sparse is dense

# lv1 is minimal sparsity
# lv2 is 0.25 ratio
# lv5 is full sparsity

p = ggplot(data=mario_time, aes(x=totalcell, y=ovlp_all, group=sparse_level)) +
  geom_line(aes(color=sparse_level))+
  geom_point(aes(color=sparse_level))+ theme_classic() + ggtitle("MARIO ovlp and all Matching Speed")
ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/sup_simulation/code_0713/data/figs6/fig2time.svg", plot=p, width=6, height=3)
## fig3 data
## first part the plot for time, fig3
mario_time = read.csv("/home/bkzhu/SNE-multi/figure_rcode/sup_simulation/code_0713/data/timeit_fig3/timeit.csv")
mario_time$totalcell = mario_time$n1 + mario_time$n2
mario_time$ovlp_all = mario_time$time_ovlp + mario_time$time_all
mario_time$sparse_level = rep(c("lv1","lv1","lv2","lv2","lv3","lv3","lv4","lv4","lv5","lv5"),5)
mario_time = subset(mario_time,mario_time$mode == "sparse") # dont need dense time since full sparse is dense
mario_time = subset(mario_time,mario_time$sparse_level %in% c("lv1","lv2","lv5") )# dont need dense time since full sparse is dense

p = ggplot(data=mario_time, aes(x=totalcell, y=ovlp_all, group=sparse_level)) +
  geom_line(aes(color=sparse_level))+
  geom_point(aes(color=sparse_level))+ theme_classic() + ggtitle("MARIO ovlp and all Matching Speed")
ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/sup_simulation/code_0713/data/figs6/fig3time.svg", plot=p, width=6, height=3)

## fig4 data
## first part the plot for time, fig4
mario_time = read.csv("/home/bkzhu/SNE-multi/figure_rcode/sup_simulation/code_0713/data/timeit_fig4/timeit.csv")
mario_time$totalcell = mario_time$n1 + mario_time$n2
mario_time$ovlp_all = mario_time$time_ovlp + mario_time$time_all
mario_time$sparse_level = rep(c("lv1","lv1","lv2","lv2","lv3","lv3","lv4","lv4","lv5","lv5"),5)
mario_time = subset(mario_time,mario_time$mode == "sparse") # dont need dense time since full sparse is dense
mario_time = subset(mario_time,mario_time$sparse_level %in% c("lv1","lv2","lv5") )# dont need dense time since full sparse is dense

p = ggplot(data=mario_time, aes(x=totalcell, y=ovlp_all, group=sparse_level)) +
  geom_line(aes(color=sparse_level))+
  geom_point(aes(color=sparse_level))+ theme_classic() + ggtitle("MARIO ovlp and all Matching Speed")
ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/sup_simulation/code_0713/data/figs6/fig4time.svg", plot=p, width=6, height=3)

############################ sparse vs matching performance ##################################

# fig 2 data
mario_sparse = read.csv("/home/bkzhu/SNE-multi/figure_rcode/sup_simulation/code_0713/data/sparsify_fig2/metrics.csv")
p = ggplot(data=mario_sparse, aes(x=sparsity, y=final_acc)) +
  geom_line()+
  geom_point()+ theme_classic() + ylim(c(0.93,0.96)) +
  ggtitle("MARIO sparsity and matching accuracy")
ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/sup_simulation/code_0713/data/figs6/fig2sparse.svg", plot=p, width=6, height=3)

# fig 3 data
mario_sparse = read.csv("/home/bkzhu/SNE-multi/figure_rcode/sup_simulation/code_0713/data/sparsify_fig3/metrics.csv")
p = ggplot(data=mario_sparse, aes(x=sparsity, y=final_acc)) +
  geom_line()+
  geom_point()+ theme_classic() + ylim(c(0.93,0.96)) +
  ggtitle("MARIO sparsity and matching accuracy")
ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/sup_simulation/code_0713/data/figs6/fig3sparse.svg", plot=p, width=6, height=3)

# fig 4 data
mario_sparse = read.csv("/home/bkzhu/SNE-multi/figure_rcode/sup_simulation/code_0713/data/sparsify_fig4/metrics.csv")
p = ggplot(data=mario_sparse, aes(x=sparsity, y=final_acc)) +
  geom_line()+
  geom_point()+ theme_classic() + ylim(c(0.85,0.86)) +
  ggtitle("MARIO sparsity and matching accuracy")
ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/sup_simulation/code_0713/data/figs6/fig4sparse.svg", plot=p, width=6, height=3)

############################ memory usage and total cell ##################################

mario_memory = readxl::read_excel("/home/bkzhu/SNE-multi/figure_rcode/sup_simulation/code_0713/data/memory_record.xlsx")
mario_memory$memory = mario_memory$memory/1024 # transfer to gb

mario_memory_fig2  = subset(mario_memory,mario_memory$data == "fig2")
p = ggplot(data=mario_memory_fig2, aes(x=`total cell`, y=memory)) +
  geom_line()+
  geom_point()+ theme_classic() + ggtitle("MARIO memory usage") +
  ggtitle("MARIO peak memory usage")
ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/sup_simulation/code_0713/data/figs6/fig2memo.svg", plot=p, width=6, height=3)

mario_memory_fig3  = subset(mario_memory,mario_memory$data == "fig3")
p = ggplot(data=mario_memory_fig2, aes(x=`total cell`, y=memory)) +
  geom_line()+
  geom_point()+ theme_classic() + ggtitle("MARIO memory usage") +
  ggtitle("MARIO peak memory usage")
ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/sup_simulation/code_0713/data/figs6/fig3memo.svg", plot=p, width=6, height=3)

mario_memory_fig4  = subset(mario_memory,mario_memory$data == "fig4")
p = ggplot(data=mario_memory_fig2, aes(x=`total cell`, y=memory)) +
  geom_line()+
  geom_point()+ theme_classic() + ggtitle("MARIO memory usage") +
  ggtitle("MARIO peak memory usage")
ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/sup_simulation/code_0713/data/figs6/fig4memo.svg", plot=p, width=6, height=3)
