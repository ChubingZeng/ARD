#' Automatic Relevance Determination Plots
#'
#' \code{plot.ARD} returns two plots -- the marginal likelihood trace plot and the coefficient plot
#' @param ARD.result ARD object
#' @return marginal likelihood trace plot and the coefficient plot likelihood scores by Automated Relevance Determination (ARD).
#' @examples
#' plot.ARD(ARD.result)

plot.ARD <- function(ARD.result){
        df = data.frame(cbind(1:ARD.result$n_iter,ARD.result$likelihood.score))
        p.likelihood<-ggplot(data=df, aes(x=X1, y=X2,group=1)) +
                geom_line(color="blue")+
                geom_point(size = 0.2)+
                theme_classic()+
                ggtitle("Marginal Log likelihood")+
                xlab('Iterations')+
                ylab("Likelihood Score")

        df.coef = data.frame(ARD.result$coefficients)
        df.relevent = data.frame(x=ARD.result$coefficients[abs(ARD.result$coefficients)>0.1],y=length(ARD.result$coefficients)/3)
        df.relevent$z <- "Features with absolute coefficient value >0.02"
        p.coefficients<-ggplot(df.coef, aes(x=ARD.result.coefficients))+
                geom_histogram(color="lightblue3", fill="deeppink1",
                               linetype="dashed",binwidth=0.02)+
                theme_classic()+
                ggtitle("Histogram of the coefficients")+
                xlab("Value of the coefficients")+
                ylab("Frequency")+
                geom_point(data = df.relevent,aes(x=x,y=y,colour = z),size=3)+
                scale_colour_manual(values = "goldenrod1",name=NULL)+
                theme(legend.background = element_rect(colour = "black",fill="lightblue1", size = 0.2),
                      legend.position = c(0.8,0.8),
                      legend.text=element_text(size=14))

        return(list(p.likelihood,p.coefficients))
}
