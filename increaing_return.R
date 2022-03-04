#2021 학부연구생 논문 개념 구현


library(plyr); library(igraph); library(reshape2); library(dplyr); library(tidyr)
#data <- read.delim("C:/Users/aaaa8889/Desktop/R/Research Program/year_origin_hs92_4.tsv", 
#                   colClasses = "character", na.strings = "NA")

data <- read.delim("~/Dropbox (Personal)/research_after 2019/2021_increasing Return/data/year_origin_hs92_4.tsv", 
                   colClasses = "character", na.strings = "NA")

#필요한 변수를 추출합니다.
export_data <- data[, c("year","origin", "hs92", "export_val")]
#데이터 전처리(NULL값을 0으로 치환)
export_data$export_val[export_data$export_val == "NULL"] <- "0"

#숫자 데이터를 알맞게 형변환하여 분석을 진행합니다.
for (j in 3:4){
  export_data[, j] <- as.numeric(as.character(export_data[, j]))
}

#str(export_data)

#RCA를 구하기 위해서 지역별, 상품별로 합계를 구합니다. 
#전체 상품간 prosimity를 구하기 위해서 시간은 고려하지 않습니다.
export_data <- aggregate(export_data$export_val, by = export_data[c("origin", "hs92")], FUN = sum)

colnames(export_data)[3] <- "export_val"

#함수 RCA에 적용하기 위해 기초적인 전처리를 시작합니다.
#행에는 origin, 열에는 hs92가 들어가고 각 셀에는 export_val이 들어갑니다.
df_val <- dcast(export_data, origin ~ hs92, value.var = "export_val", fill = 0)
View(df_val)
#함수 적용을 위해 행렬로 변환합니다. 원하는 변수가 2번째부터 있기 때문에 2:ncol()을 사용하였습니다.
mat_val <- as.matrix(df_val[,2:ncol(df_val)])

#RCA를 계산하는 함수를 선언합니다. package "EconGeo"에서 가져왔습니다.
RCA <- function (mat, binary = FALSE) {
  mat <- as.matrix(mat)
  share_tech_city  <- mat/rowSums(mat)
  share_tech_total <- colSums(mat)/sum(mat)
  
  LQ <- t(t(share_tech_city)/share_tech_total)
  LQ[is.na(LQ)] <- 0
  
  if (binary) {
    LQ[LQ < 1] <- 0
    LQ[LQ > 1] <- 1
  } 
  LQ = round (LQ, digits = 2)
  return(LQ)
}

#변환된 행렬을 가지고 RCA를 구합니다.
mat_rca <- RCA(mat_val)


#행이름에 origin을 부여합니다.
rownames(mat_rca) <- df_val[ , 1]

#계산된 RCA를 바탕으로 1보다 크면 1을, 1보다 작으면 0을 리턴합니다.
mat_mcp <- ifelse(mat_rca > 1, 1, 0)

#전치행렬과 행렬을 곱하여 product by product 행렬을 구합니다. 
#행렬의 원소는 상품간 교집합 RCA 수를 나타냅니다.
mat_multi <- t(mat_mcp)%*%mat_mcp

#proximity를 계산하는 함수를 선언합니다.
#자세한 내용은 코드정리를 참고해주세요.
make_proximity <- function(mat){
  # mat is a product by product matrix of which element is the number of country that has RCA(i) over 1.
  # 대각원소를 선택합니다.
  diag_vector <- diag(mat)
  
  #대각원소를 정렬하여 계산을 용이하게 합니다.
  e <- order(diag_vector, decreasing = T)
  
  # 연산의 바탕이 되는 행렬을 저장합니다.
  l = mat
  
  #연산된 결과를 저장하는 행렬을 저장합니다.
  b = mat
  for (k in e){
    #대각벡터에서 k번째 해당하는 원소를 꺼내 l의 k번째 행에다 나누어줍니다.
    #그리고 그 값이 연산의 결과가 저장된 b행렬과 서로 비교하여 벡터의 원소끼리 대소비교를 합니다.
    p <- (b[k, ] > (1/diag_vector[k])*l[k, ])
    
    #대소비교를 해서 True값이 나온 값에다가 대각원소로 나눈 원소를 할당합니다.
    #이 과정을 통해 Proximity의 조건부 확률을 잡아냅니다.
    b[k, p] <- (1/diag_vector[k])*l[k, p]
    
    #위의 연산을 열에다가 반복하는 작업입니다.
    q <- b[,k] > (1/diag_vector[k])*l[, k]
    #마찬가지로 위 연산을 열에데가 반복하는 작업입니다.
    b[q,k] <- (1/diag_vector[k])*l[q, k]
  }
  return(b)
}

#product by product 행렬을 가지고 product space를 생성합니다.
proximity <- make_proximity(mat_multi)


#proximity를 기반으로 각 국가 및 상품의 dencity를 구합니다.
#dencity의 경우 시간에 따라 변하는 특성을 가지므로 시간축을 따로 설정합니다.
#unique()함수를 사용하여 시간축을 설정합니다.
year_set <- unique(data$year)

#시간축을 충심으로 다음의 계산이 진행됩니다.
for (i in year_set){
  #필요한 변수를 추출합니다. 
  #line 42번 코드와 다른점은 이번 변수추출은 각 국가의 시간별 RCA를 구하기 위해서 입니다.
  export_data <- data[data$year == i, c("year","origin", "hs92", "export_rca")]
  export_data$export_rca[export_data$export_rca == "NULL"] <- "0"
  for (j in 3:4){
    export_data[, j] <- as.numeric(as.character(export_data[, j]))
  }
  #Mcp 변수를 선언합니다.
  export_data$Mcp <- 0
  #export_rca가 1보다 큰 경우 1을 부여합니다. 그러면 Mcp변수는 자연스럽게 0,1를 가집니다.
  export_data[which(export_data$export_rca >= 1), 5] <- 1
  #Mcp가 1인 데이터프레임을 추출합니다.
  ex <- export_data[which(export_data$Mcp == 1),]
  
  
  m <- ex[, c("origin", "hs92")]
  m[,2] <- as.numeric(as.character(m[,2]))
  
  mydata <- melt(m, id.vars = "origin")
  #View(mydata)
  M <- dcast(mydata, origin ~ value, function(x) 1, fill = 0)
  #View(M)
  Mcp <- as.matrix(M[,2:ncol(M)])
  #View(Mcp)
  rownames(Mcp)<- M[,1]
  #Mcp_multip <- t(Mcp)%*%Mcp
  
  
  #Mcp 행렬과 proximity의 행렬의 차원이 서로 다르기 때문에 차원을 일치시킵니다.
  #Mcp에 proximity의 상품명이 없는지 판단하는 코드입니다.
  #코드 결과로 k에는 proximity에는 있으나, mcp에는 없는 상품명이 추출됩니다.
  k <- !(colnames(proximity) %in% colnames(Mcp))
  
  #상품명을 추가하여 차원을 동일하게 맞춥니다.
  for (h in colnames(proximity)[k]){
    #우선 새로운 0 벡터를 mcp행렬에 추가하고
    Mcp <-cbind(Mcp, 0)
    #새롭게 생성된 마지막 열이름을 h로 설정합니다.
    colnames(Mcp)[-(1:length(colnames(Mcp))-1)] <- h
  }
  
  #추가된 열이름을 proximity 열 이름과 동일하게 정렬합니다.
  Mcp <- Mcp[, order(as.numeric(colnames(Mcp)))]
  
  #Density를 계산합니다.
  #Density의 정의에 따라 상품마다 proximity합을 구합니다.
  degrees <- colSums(proximity)
  
  #이를 데이터 프레임으로 변환합니다.
  degrees <- data.frame(degrees)
  
  #degrees에서 상품명을 추출합니다.
  degrees$product <- rownames(degrees)
  
  #density 정의대로 계산해줍니다.
  densityDeveloped = Mcp%*%proximity
  densityAroundProducts  <- densityDeveloped/t(replicate(dim(densityDeveloped)[1],rowSums (proximity,dims = 1)))
  
  #정의식대로 구한 결과를 데이터프레임으로 변환합니다.
  densityAroundProducts2 <- as.data.frame(densityAroundProducts)
  densityAroundProducts2$origin <- rownames(densityAroundProducts2)
  
  #merge를 위해 origin을 기준으로 데이터를 변형합니다.
  molten <- melt(densityAroundProducts2, id.vars = "origin")
  colnames(molten) <- c("origin", "hs92", "omega")
  molten$hs92 <- as.numeric(as.character(molten$hs92))
  
  #위 계산으로 얻어진 omega를 origin, hs92를 기준으로 병합합니다.
  trade_data <- left_join(export_data, molten, by = c("origin", "hs92"))
  
  #기본적으로 1995년도 데이터를 생성한다음
  if (i == "1995"){
    set_trade_data <- trade_data
    
    #이후 년도부터는 계속해서 rbind하여 데이터를 이어붙입니다.  
  } else{
    set_trade_data <- rbind(set_trade_data, trade_data)
  }
}

print("end")






#계산편의를 위해 다시 할당합니다.

refine <- set_trade_data
data$export_val[data$export_val == "NULL"] <- 0
data$export_val <- as.numeric(data$export_val)
#시기별 국가별 trade 합을 구합니다.
trade_sum <- aggregate(data$export_val, by = data[c("year", "origin")], FUN = sum)
colnames(trade_sum) <- c("year", "origin", "trade_val")

#시간축과 국가별로 데이터를 병합하여 trade총계를 병합합니다.
nr <- merge(refine, trade_sum, by = c("year", "origin"))

###trade_t와 trade_t+2변수를 생성하는 과정입니다.
nr$year <- as.numeric(nr$year)
for (i in 1995:2015){
  #t기 trade 데이터를 추출하고
  t_data <- nr[nr$year == i, ]
  #t+2기 trade 데이터를 추출합니다.
  next_data <- nr[nr$year == i+2, c("origin", "hs92", "Mcp", "trade_val")]
  #두개의 데이터프레임을 origin, hs92로 병합합니다.
  t_data <- merge(t_data, next_data, by = c("origin", "hs92"))
  
  if(i == 1995){
    r <- t_data
  } else{
    r <- rbind(r, t_data)
  }
  
  
}

print("end")

#보기 좋게 이름을 다시 설정합니다.
colnames(r)[5] <- "Mcp_t"
colnames(r)[7:9] <- c("trade_t", "Mcp_t2", "trade_t2")

final_data <- r
final_data$year <- as.character(final_data$year)
###

#데이터를 불러옵니다.
gdp <- read.csv("GDP.csv", colClasses = "character", na.strings = "NA")
capital <- read.csv("capital.csv", colClasses = "character", na.strings = "NA")
pop <- read.csv("pop.csv", colClasses = "character", na.strings = "NA")
saving <- read.csv("saving.csv", colClasses = "character", na.strings = "NA")
for_inflow <- read.csv("for_inflow.csv", colClasses = "character", na.strings = "NA")
export_size <- read.csv("export_size.csv", colClasses = "character", na.strings = "NA")

#데이터를 합칠 수 있도록 데이터구조를 변형합니다.
make_data <- function(x){
  x <- x[, 2:23]
  colname <- 1995:2015
  colname <- as.character(colname)
  colnames(x)[2:22] <- colname
  #View(x)
  
  x_melt <- melt(x, id.vars = "Country.Code")
  colnames(x_melt) <- c("origin", "year", "export_size")
  x <- x_melt
  
  x$origin <- tolower(x$origin)
  x$year <- as.character(x$year)
  
  return(x)
}
gdp <- make_data(gdp)
capital <- make_data(capital)
pop <- make_data(pop)
saving <- make_data(saving)
for_inflow <- make_data(for_inflow)
export_size <- make_data(export_size)

#데이터를 origin과 year로 병합합니다.
final_data <- left_join(final_data, gdp, by = c("origin", "year"))
final_data <- left_join(final_data, saving, by = c("origin", "year"))
final_data <- left_join(final_data, pop, by = c("origin", "year"))
final_data <- left_join(final_data, capital, by = c("origin", "year"))
final_data <- left_join(final_data, for_inflow, by = c("origin", "year"))

#완성된 데이터 입니다.
View(final_data)

final_data <- final_data[c(3,5,8,7,9,1,2,4,6,10,11,12,13,14)]

assignment <- read.csv("assignment.csv", colClasses = "character",na.strings="NA")

for (i in c(2:6, 8:15)) {
  assignment[,i] <- as.numeric(as.character(assignment[,i]))}

assignment$square <- assignment$omega * assignment$omega

trade <- assignment[which(assignment$Mcp_t != 1), ]

trade <- trade[which(trade$year >= 2012),]

stat1 <- lm(trade_t2 ~ log(omega+0.0001) + trade_t ,data= trade)
stat2 <- lm(trade_t2 ~ log(omega+0.0001) + trade_t + log(square + 0.0001), data= trade)
stat3 <- lm(trade_t2 ~ log(omega+0.0001) + trade_t + log(square + 0.0001) + gdp + pop, data= trade)
summary(stat3)

stat1 <- glm(Mcp_t2  ~ log(omega+0.0001)  , family = "binomial", data= trade)
stat2 <- glm(Mcp_t2  ~ log(omega+0.0001)  + log(square + 0.0001), family = "binomial", data= trade)
stat3 <- glm(Mcp_t2  ~ log(omega+0.0001)  + log(square + 0.0001) + gdp + pop, family = "binomial",data= trade)
summary(stat3)

stargazer(stat1, stat2, stat3, out = 'trade.html' )

###End