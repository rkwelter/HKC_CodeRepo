function [dX] =HKC1_LaTeX(X,Pr,Ra,Ro,k1,V)

\begin{equation} X = \begin{pmatrix} u^+_{(0,1)} \\
u^+_{(1,1)} \\
u^-_{(0,1)} \\
u^-_{(1,1)} \\
\theta_{(0,2)} \\
\theta_{(1,1)} \\
\end{pmatrix} \end{equation}

\begin{equation} \begin{split} \frac{d}{dt} u^+_{(0,1)} & =-\Pra u^+_{(0,1)} + \Pra \Rot u^-_{(0,1)} \\
\frac{d}{dt} u^+_{(1,1)} & =-(k_1^2 +1)\Pra u^+_{(1,1)} - \frac{1\Pra \Ray k_1}{\sqrt{k_1^2 +1}} \theta_{(1,1)} + \frac{1\Pra \Rot }{\sqrt{k_1^2 +1}} u^-_{(1,1)} \\
\frac{d}{dt} u^-_{(0,1)} & =-\Pra u^-_{(0,1)} - \Pra \Rot u^+_{(0,1)} \\
\frac{d}{dt} u^-_{(1,1)} & =-(k_1^2 +1)\Pra u^-_{(1,1)} - \frac{1\Pra \Rot }{\sqrt{k_1^2 +1}} u^-_{(1,1)} \\
\frac{d}{dt} \theta_{(0,2)} & =-4 \theta_{(0,2)} - \frac{k_1}{4V} \Big [ \frac{2^{\frac{1}{2}}(-2)}{\sqrt{k_1^2 +1}} u^+_{(1,1)} \theta_{(1,1)} \Big ] \\
\frac{d}{dt} \theta_{(1,1)} & =-(k_1^2 +1) \theta_{(1,1)} - \frac{1 k_1}{\sqrt{k_1^2 +1}} u^+_{(1,1)} - \frac{k_1}{4V} \Big [ \frac{2^{\frac{1}{2}}(2)}{\sqrt{k_1^2 +1}} u^+_{(1,1)} \theta_{(0,2)} \Big ]\end{split} \end{equation}
end