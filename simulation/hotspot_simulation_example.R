#' ---
#' title: "HotSpot Simulation and Analysis"
#' author: "Jiaqiang Zhu"
#' date: "Oct 15th, 2020"
#' ---


#' >  Generate the data  
#+ label=dg, echo=T, warnings=F, message=F,eval=F

rm(list=ls())
library(spatstat)
library(Matrix)
library(MASS)

source("./funcs/utilities.R")

numGenes = 10000 # 500
lowSignal = 500 # 100
highSignal = 500

numCell = 3000
cellProp = 0.2## 0.1, 0.2, 0.3
irpt = 1
mu0 = 0.5
ieffect = 3
itheta = 5

cat("\n ## of Cells",numCell,"\n")
cat(irpt," ")
set.seed(irpt)
pp 			<- rpoispp(numCell)
low_expr 	<- c(1)
high_expr 	<- c(2)

pp = add_markdist_hotspot(pp,low_expr, high_expr,cell_proportion=cellProp)
back_ground_cells 	<- which(pp$marks==1)
spiked_in_cells 	<- which(pp$marks==2)


expr_mat <- matrix(NA,nrow=pp$n,ncol=numGenes)
for(igene in 1:numGenes){
	a <- rnegbin(pp$n,mu0,theta=itheta)
	if(igene <= lowSignal){
		a[spiked_in_cells] 	<- rnegbin(length(spiked_in_cells),mu0*ieffect,theta=itheta)
	}else if(igene <=lowSignal+highSignal) {
		a[spiked_in_cells]  <- rnegbin(length(spiked_in_cells),mu0/ieffect,theta=itheta)
	}
	expr_mat[,igene] <- a
}

colnames(expr_mat) 	<- paste0("gene",1:numGenes)
info 				<- cbind.data.frame(x=pp$x,y=pp$y)
rownames(info) 		<- rownames(expr_mat) <- paste0("C",1:nrow(expr_mat))

sp_sim_count 		<- as(t(expr_mat),"sparseMatrix")
location 			<- info 

save(sp_sim_count,location,file=paste0("./simulation/sim_hotspot_SS",numCell,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_rpt",irpt,".rds"))

		

#' *** 
#' >  Check Gene Patterns 
#+ label=pattern,fig.width=12, fig.height=4,warnings=F, echo=T,fig.align="center",eval=T
rm(list=ls())
source("./funcs/utilities.R")
numCell = 3000
cellProp = 0.2
irpt = 1
mu0 = 0.5
ieffect = 3
itheta = 5
load(paste0("./simulation/sim_hotspot_SS",numCell,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_rpt",irpt,".rds"))
pltdat  <- cbind.data.frame(location[,1:2],apply(sp_sim_count[c(1,501,1001),],1,relative_func))
pp 		<- lapply(1:3,function(x){pattern_plot_sparkx(pltdat,x,main=F,pointsize=1.5,titlesize=1.5,
								min.pand=0.8,max.pand=1.01,title="Hotspot", 
								pal=colorRampPalette(c("antiquewhite",
									viridis_pal()(10)[c(6,2)])))})
library(ggpubr)
fig1 	<- ggarrange(pp[[1]], pp[[2]],pp[[3]],ncol = 3, nrow = 1)
fig1




#' *** 
#' >  Analyze the data with SPARK-X
#+ label=analysis,warnings=F,echo=T,eval=T
rm(list=ls())
library(SPARK)
numCell = 3000
cellProp = 0.2
irpt = 1
mu0 = 0.5
ieffect = 3
itheta = 5
load(paste0("./simulation/sim_hotspot_SS",numCell,"_cellProp",cellProp,"_mu",mu0,"_effect",ieffect,"_theta",itheta,"_rpt",irpt,".rds"))
sparkX <- sparkx(sp_sim_count,as.matrix(location))
head(sparkX$res_mtest)



