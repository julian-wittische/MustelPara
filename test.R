library(performance)
library(glmmTMB)
data(Owls)

m_RAND <- glmmTMB(NegPerChick~BroodSize + ArrivalTime + (1|Nest), data=Owls)
m_NORAND <- glmmTMB(NegPerChick~BroodSize + ArrivalTime, data=Owls)

r2(m_RAND)
r.squaredGLMM(m_RAND)

r2(m_NORAND)
r.squaredGLMM(m_NORAND)


