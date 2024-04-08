
## This script shows some examples of plotting with ggplot2, base R and plotly.
#Can be used to either transform ggplot2 to interactive plots or to directly create interactive plots using plotly functions

# Load Libraries
library(ggplot2)
library(dplyr)
library(plotly)
library(hrbrthemes)
library(htmlwidgets)
library(tidyverse)
library(viridis)
library(heatmaply)

#Load sample data (iris dataset)
data(iris)
summary(iris)
dim(iris)
variables <- names(iris)


#### Scatterplot

## ggplot2

ggplot.1 <- ggplot(iris, aes(x=Petal.Length, y=Petal.Width)) +
  geom_point(color="#FFA500", alpha=0.8, size=3) +
  ggtitle("Petal length and width") +
  #theme_ipsum() + optional theme
  theme(plot.title = element_text(size=15)) +
  ylab('Petal width') +
  xlab('Petal length')
ggplot.1

## Base R

plot(iris$Petal.Length, iris$Petal.Width, main = "Petal length and width",
     xlab = 'Petal length', ylab = 'Petal width',
     pch = 19, col=alpha("#FFA500", 0.8))
grid(nx = 11, ny = 8)

## Plotly

plotly.1 <- plot_ly(data = iris, x = ~Petal.Length, y = ~Petal.Width,type="scatter", opacity=0.7,mode="markers",text = ~Species,size = 3,
                    marker=list(color="#FFA500",line = list(color = "#FFA500",
                                                             width = 2)))%>% 
  layout(title = "Petal length and width",
         yaxis = list(title='Petal Length'),
         xaxis = list(title='Petal width'))
plotly.1

#### Multi-histogram

##ggplot2
ggplot.2 <- ggplot(iris, aes(x=Sepal.Width, fill=Species)) +
  geom_histogram(color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#FFA500","#404080" )) +
  #theme_ipsum() +
  theme(plot.title = element_text(size=15)) +
  ggtitle("Sepal width")+
  ylab('Count') +
  xlab('Sepal width')
ggplot.2

##base R

setosa <- iris[iris$Species=="setosa",]$Sepal.Width
versicol <- iris[iris$Species=="versicolor",]$Sepal.Width
virginica <- iris[iris$Species=="virginica",]$Sepal.Width

grid()
# First distribution
hist(setosa, breaks=30,border="#e9ecef", col=alpha("#69b3a2",0),xlab="height", 
     ylab="Count", main="Histogram" )
grid(nx = 11, ny = 8)
hist(setosa, breaks=30, border="#e9ecef",col=alpha("#69b3a2", 0.6), xlab="height", 
     ylab="Count", main="Histogram",add=T)
# Second with add=T to plot on top
hist(versicol, breaks=30, col=alpha("#404080",0.6),border="#e9ecef",add=T)
# Third with add=T to plot on top
hist(virginica, breaks=30, col=alpha("#FFA500",0.6),border="#e9ecef",add=T)
# Add legend
legend("topright", legend=c("Setosa","Versicolor", "Virginica"), col=c("#69b3a2","#404080","#FFA500"), pt.cex=2, pch=15 )


## plotly 
plotly.2 <- plot_ly(alpha = 0.6) %>% 
  add_histogram(x = setosa, opacity=0.6,marker= list(color ="#404080",line =list(color="#e9ecef", width=1)), name= "Setosa") %>% 
  add_histogram(x = versicol,opacity=0.6,marker= list(color = "#69b3a2",line =list(color="#e9ecef", width=1)),name= "Versicolor") %>% 
  add_histogram(x = virginica,opacity=0.6,marker= list(color = "#FFA500",line =list(color="#e9ecef", width=1)),name= "Virginica") %>% 
  
  layout(barmode = "overlay", title="Histogram example", xaxis = list(title = "Value",showgrid = TRUE), yaxis = list(title = "Count",showgrid = TRUE)
  )
plotly.2 

##---------------------------

#Example 1: Area chart
p <- data %>%
  ggplot( aes(x=date, y=value)) +
  geom_area(fill="#69b3a2", alpha=0.5) +
  geom_line(color="#69b3a2") +
  ylab("bitcoin price ($)") +
  theme_ipsum()

# Turn it interactive with ggplotly
p <- ggplotly(p)
p

# save the widget: it is saved as a web link to interact with the plot created.
saveWidget(p, file=paste0(getwd(), "/ggplotlyAreachart.html"))

#Example 2: 
# Create data
var1 <- seq(1,8)
var2 <- c(0,1,4,1,8,7,5,4)
var3 <- c(7,8,4,2,1,2,0,1)

# Area chart with 2 groups
p <- plot_ly(x = var1, y = var2, type="scatter", mode="markers", fill = "tozeroy")
p <- add_trace(p, x = var1, y = var3, type="scatter", mode="markers", fill = "tonexty")
p

# save the widget
# saveWidget(p, file=paste0( getwd(), "/HtmlWidget/plotlyAreachart.html"))

#Example 3: Interactive heatmap

# Load data 
data <- read.table("https://raw.githubusercontent.com/holtzy/data_to_viz/master/Example_dataset/multivariate.csv", header=T, sep=";")
colnames(data) <- gsub("\\.", " ", colnames(data))

# Select a few country
data <- data %>% 
  filter(Country %in% c("France", "Sweden", "Italy", "Spain", "England", "Portugal", "Greece", "Peru", "Chile", "Brazil", "Argentina", "Bolivia", "Venezuela", "Australia", "New Zealand", "Fiji", "China", "India", "Thailand", "Afghanistan", "Bangladesh", "United States of America", "Canada", "Burundi", "Angola", "Kenya", "Togo")) %>%
  arrange(Country) %>%
  mutate(Country = factor(Country, Country))

# Matrix format
mat <- data
rownames(mat) <- mat[,1]
mat <- mat %>% dplyr::select(-Country, -Group, -Continent)
mat <- as.matrix(mat)


p <- heatmaply(mat, 
               dendrogram = "none",
               xlab = "", ylab = "", 
               main = "",
               scale = "column",
               margins = c(60,100,40,20),
               grid_color = "white",
               grid_width = 0.00001,
               titleX = FALSE,
               hide_colorbar = TRUE,
               branches_lwd = 0.1,
               label_names = c("Country", "Feature:", "Value"),
               fontsize_row = 5, fontsize_col = 5,
               labCol = colnames(mat),
               labRow = rownames(mat),
               heatmap_layers = theme(axis.line=element_blank())
)

#Example 4: Scatterplot

data <- read.table("https://raw.githubusercontent.com/holtzy/data_to_viz/master/Example_dataset/2_TwoNum.csv", header=T, sep=",") %>% dplyr::select(GrLivArea, SalePrice)

##ggplot2
plot1 <- ggplot( data, aes(x=GrLivArea, y=SalePrice/1000)) +
  geom_point(color="#69b3a2", alpha=0.8) +
  ggtitle("Ground living area partially explains sale price of apartments") +
  theme_ipsum() +
  theme(plot.title = element_text(size=12)) +
  ylab('Sale price (k$)') +
  xlab('Ground living area')
plot1
##base R

plot(data$GrLivArea, data$SalePrice/1000, main = "Ground living area partially explains sale price of apartments",
     xlab = 'Ground living area', ylab = 'Sale price (k$)',
     pch = 19, col=alpha("#69b3a2", 0.8),frame = FALSE)
grid(nx = 11, ny = 8)

## Plotly

fig <- plot_ly(data = data, x = ~GrLivArea, y = ~SalePrice/1000,type="scatter", opacity=0.7,mode="markers", marker=list(color="#69b3a2"))%>% 
  layout(title = "Ground living area partially explains sale price of apartments",
           yaxis = list(title='Sale price (k$)'),
           xaxis = list(title='Ground living area'),showgrid = TRUE)
fig

#Example 5: Multi-Histogram

data <- data.frame(
  type = c( rep("variable 1", 1000), rep("variable 2", 1000) ),
  value = c( rnorm(1000), rnorm(1000, mean=4) )
)

##ggplot2
plot2 <- ggplot(data,aes(x=value, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  theme_ipsum() +
  labs(fill="")
plot2

## base R

h1 <- data[data$type=="variable 1",]$value
h2 <- data[data$type=="variable 2",]$value

grid()
# First distribution
hist(h1, breaks=30, xlim=c(min(data$value),max(data$value)),border="#e9ecef", col=alpha("#69b3a2",0),xlab="height", 
     ylab="Count", main="Histogram" )
grid(nx = 11, ny = 8)
hist(h1, breaks=30, xlim=c(min(data$value),max(data$value)),border="#e9ecef",col=alpha("#404080", 0.6), xlab="height", 
     ylab="Count", main="Histogram",add=T)
# Second with add=T to plot on top
hist(h2, breaks=30, col=alpha("#69b3a2",0.6),border="#e9ecef",add=T)

# Add legend
legend("topright", legend=c("Variable 1","Variable 2"), col=c("#404080", 
                                                              "#69b3a2"), pt.cex=2, pch=15 )

#plotly 
fig <- plot_ly(alpha = 0.6) %>% 
  add_histogram(x = h1, opacity=0.6,marker= list(color ="#404080",line =list(color="#e9ecef", width=1)), name= "Variable 1") %>% 
  add_histogram(x = h2,opacity=0.6,marker= list(color = "#69b3a2",line =list(color="#e9ecef", width=1)),name= "Variable 2") %>% 
  layout(barmode = "overlay", title="Histogram example", xaxis = list(title = "Value",showgrid = TRUE), yaxis = list(title = "Count",showgrid = TRUE)
  )
fig  

## Example 6: Grouped boxplot

variety=rep(LETTERS[1:7], each=40)
treatment=rep(c("high","low"),each=20)
note=seq(1:280)+sample(1:150, 280, replace=T)
data=data.frame(variety, treatment ,  note)

# ggplot2

plot2 <- ggplot(data, aes(x=variety, y=note, fill=treatment)) + 
  geom_boxplot()+ labs(title="Grouped boxplot with ggplot2")+scale_fill_manual(values=c("slateblue1" , "tomato"))
plot2

# Base R
par(new = TRUE)
bop<-boxplot(note~treatment*variety, data=data,
        col = c("slateblue1" , "tomato"), 
        main = "Grouped Boxplot with base R", 
        xlab = "Variety", 
        ylab = "Note",   
        outline = TRUE,    # Optional: Show outliers
        horizontal = FALSE,xaxt = "n")
my_names <- sapply(strsplit(bop$names , '\\.') , function(x) x[[2]] )
my_names <- my_names[seq(1 , length(my_names) , 2)]
bop <- axis(1, 
     at = seq(1.5 , 14 , 2), 
     labels = my_names, 
     tick=TRUE , cex=0.3)
bop<-legend("topright", legend = c("High treatment", "Low treatment"), 
       col=c("slateblue1" , "tomato"),
       pch = 15, bty = "n", pt.cex = 3, cex = 1.2,  horiz = F, inset = c(0.1, 0.1))
bop
