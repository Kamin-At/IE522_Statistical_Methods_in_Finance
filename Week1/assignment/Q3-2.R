#install.packages("ggplot2")
#install.packages("ggrepel")
library(ggplot2)
library(ggrepel)

data = c(1, 2, 4, 6, 6.75, 9, 12)
data_size = length(data) - 2
Quantile = c(0/data_size, 
          1/data_size, 
          2/data_size, 
          3/data_size, 
          0.65, 
          4/data_size, 
          5/data_size)

df = data.frame(x=data,
                y=Quantile)
z = Map( function(x, y){toString(c(x, y))}, df$x, df$y)
print(paste('Z:', z))
print(df)

ggplot(data=df, aes(x=x, y=y, group=1)) +
geom_line()+
geom_point()+
geom_text_repel(aes(label = z))