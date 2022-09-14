##  Análise de Resíduos_________________________________________________________
##  Objetivo: verificar a adequação da distribuição escolhida e identificar possiveis outliers.

#   Ajuste dos Modelo Lineares Generalizados____________________________________
#   Dados simulados

set.seed(1)

## amostra de tamanho 30
n = 30

## Variaveis explicativas que recebem valor 0 , 1 , 2

x1 <- rep((0:2), n/3)
x2 <- rep((0:2), each = n/3)

# Valor esperado variavel y

lambda = round(exp(1.8 + 0.5 * x1 + 0.2 * x2),2)
head(lambda)

#Variavel resposta gerada a partir de uma distribuição de poisson com valor esperado

resposta = rpois(n, lambda)
dados <- data.frame(resposta, x1, x2)
dados

# Ajustando o modelo linear generealizado_______________________________________

modelo <- glm(resposta ~ x1 + x2, data = dados, 
              family = poisson(link = "log"))

summary(modelo)

## Resp:  
## intercepto bem proximo do simulado que foi 1.8
## x1 (beta0) bem proximo do simulado que foi 0.5
## x2 (beta1) bem proximo do simulado que foi 0.2
## 27 grau de liberade  (n - p = 30 - 3 parametros = 27)
# Residual deviance do modelo

## Medidas de Ajuste do modelo__________________________________________________
# Deviance = Soma das logverossimilhanas de todas observações

deviance = with(modelo, cbind(Deviance = deviance, 
                              "Graus de Liberdade" = df.residual,
                              "P-valor" = pchisq(deviance, df.residual, 
                                                 lower.tail=FALSE)))
deviance
## H0 : O modelo esta bem ajustado
## H1: O modelo não está bem ajustado

##R: p-valor > 0,05, não rejeitamos a hipotese que o modelo está bem ajustado
## ou
modelo$deviance
## graus de liberdade
modelo$df.residual

## p-valor da deviance considerando a dist qui-quadrado
1 - pchisq(modelo$deviance, modelo$df.residual)

# AIC do modelo
modelo$aic


## Matriz H_____________________________________________________________________

retorna.H <- function(modelo){
  # extrai a matriz x
  X = model.matrix(modelo)
  # extrai a diagonal da matriz w
  W = diag(modelo$weights)
  M = solve(t(X)%*%W%*%X)
  H = sqrt(W)%*%X%*%M%*%t(X)%*%sqrt(W)
  
  # Extrai os elementos da diagonal da Matriz H
  
  h = diag(H)
  return(h)}

## Matriz H

h = retorna.H(modelo)
phi = summary(modelo)$dispersion

## Desvios Residuais

r.desvio = residuals(modelo, type = "deviance")
head(r.desvio)

## Calculando de maneira manual

sinal.func <-function(y, y_hat){
  sinal = rep(1, length(y))
  sinal[y < y_hat] = -1
  return(sinal)}

y = dados$resposta
y_hat = fitted(modelo)

sinal = sinal.func(y, y_hat)
r.desvio.manual = sinal * sqrt(2 * (y * log(y/y_hat) - y + y_hat))

## Comparando os resultados

head(r.desvio.manual)
head(r.desvio)

## Calculando o Desvio Residula Padronizado

r.desvio.p = r.desvio * sqrt(phi/(1-h))


# Resíduos de Pearson___________________________________________________________

r.pearson  = residuals(modelo, type = "pearson")
head(r.pearson)

## Calculando de maneira manual residuos de Pearson

head(dados)
head(fitted(modelo))
r.pearson.manual = (dados$resposta - fitted(modelo))/sqrt(fitted(modelo))

## Comparando os resultados

head(r.pearson.manual)
head(r.pearson)

## Calculando o Resíduo de Pearson Padronizado
r.pearson.p = r.pearson * sqrt(phi/(1-h))

## Análise de Resíduos__________________________________________________________
## Tem como objetivo verificar a adequação da distribuição escolhida e identificar possiveis outliers

par(mfrow=c(2,2))
phi = summary(modelo)$dispersion
res <- residuals(modelo, type = "deviance") * sqrt(phi/(1-h));

## Gráfico de Residuos X Valores Ajustados

plot(modelo$fitted, res, xlab = "Valores ajustados",
     ylab = "Desvio Residual padronizado", pch = 19, col = "black")
abline(h = 0, col = 'red', lty = 2)
##  Obs: Em principio prece existir uma nuvem aleatória em torno de 0, desta maneira consideramos ok.

# Gráfico de Residuos X Ordem
plot(res, ylab = "Desvio Residual padronizado", xlab = "Ordem",
     pch = 19, col = "black"); abline(h = 0, col = 'red', lty = 2)
##  Obs: Em principio prece existir uma nuvem aleatória em torno de 0, desta maneira consideramos ok.


# Gráfico histograma dos resíduos
# obs: Analisar a distribuição simetrica dos residuos, aparentemente não simetrico
hist(res, main = "", ylab = "Frequência", 
     xlab = "Desvio Residual padronizado")

# Gráfico Quantil da Distribuição Normal
qqnorm(res)

# Teste de normalidado dos residuos
shapiro.test(res)	

## Simulação Monte Carlo________________________________________________________
## obs: Vereificação  da adequação do modelo a partir dos residuos e construção 
## do gráfico de envelope.

phi = summary(modelo)$dispersion
residuos.modelo = resid(modelo, type ="deviance")
residuos.modelo = residuos.modelo * sqrt(phi/(1-h))
n <- nrow(dados)
m = 100
matriz.residuos <- matrix(NA,n,m)

for(i in 1:m){
  nresp <- rpois(n, fitted(modelo))
  fit <- glm(nresp ~ x1 + x2, family = poisson(link = "log"))
  h = retorna.H(fit)
  phi = summary(modelo)$dispersion
  matriz.residuos[,i] = 
    sort(resid(fit, type = "deviance") * sqrt(phi/(1-h)))}

head(matriz.residuos)

## Cálculo dos percentis
intervalos=apply(matriz.residuos, 1, function(x) 
  quantile(x, c(0.025, .975)))

## Cálculo da média
med <- apply(matriz.residuos, 1, mean)

## Amplitude de variação dos resíduos
faixa = range(residuos.modelo, intervalos)


## Gráfico dos quantis téoridos da distribução normal
par(mfrow = c(1,1))
qqnorm(residuos.modelo, xlab = "Percentil da Normal Padrão",
       ylab = "Desvio Residual Padronizado", ylim = faixa, pch = 16)

## Adicionando o envelope simulado ao gráfico
par(new = T)
qqnorm(intervalos[1,], axes = F, xlab = "", ylab = "",type = "l",
       ylim = faixa, lty = 1, cex = 0.6)
par(new = T)
qqnorm(intervalos[2,], axes = F, xlab = "", ylab = "", type = "l",
       ylim = faixa, lty = 1, cex = 0.6)
par(new = T)
qqnorm(med, axes = F, xlab = "", ylab = "", type = "l",
       ylim = faixa, lty = 2, cex = 0.6)

