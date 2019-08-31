import msvcrt
import numpy as np
import math
FMAD = 22.7
TEAD = 33.9
TSAD = 29.4
DAD = 995
CPAD = 4.18
ViAD = 0.00082
KAD = 0.62
FMA = 35.32
TEA = 23.9
TSA = 26.7
DA = 999
CPA = 4.18
ViA = 0.00092
KA = 0.62
tm = (TEA + TSA) / 2
# Para el agua destilada
PrAD = CPAD * ViAD * 1000 / KAD
# Para el agua
PrA = CPA * ViA * 1000 / KA
t1 = TEA
t2 = TSA
T1 = TEAD
T2 = TSAD
deltaTlm = ((T1 - t2) - (T2 - t1)) / math.log((T1 - t2) / (T2 - t1))
R = (T1 - T2) / (t2 - t1)
S = (t2 - t1) / (T1 - t1)
a = ((R**2) + 1)**(1 / 2)
b = math.log((1 - S) / (1 - R * S))
c = R - 1
d = 2 - S * (R + 1 - a)
e = 2 - S * (R + 1 + a)
ft = a * b / (c * math.log(d / e))
deltaTm = ft * deltaTlm
Q = (FMAD) * (CPAD) * (T1 - T2) * 1000
ListaNp = [1, 2, 4, 6, 8]
ListaL = [1.83, 2.44, 3.66, 4.88, 6.10, 7.32]
ListaDE = [0.016, 0.020, 0.025, 0.030, 0.038, 0.050]
ListaG = [0.0016, 0.0012, 0.0020, 0.0026, 0.0032]
ListaCB = [0.15, 0.25, 0.35, 0.45]
ListaArreglo = ['Triangular', 'Cuadrado']
ListaClaroCoraza = ['Pull-through floating head', 'Split-ring floating head', 'Outside packed head', 'Fixen and U-tube']
ListaU0 = [np.arange(1500, 1501, 1)]
i = 1
for U0 in ListaU0:
	for Np in ListaNp:
		for L in ListaL:
		    for DE in ListaDE:
			    for G in ListaG:
				    for CB in ListaCB:
					    for Arreglo in ListaArreglo:
						    for ClaroCoraza in ListaClaroCoraza:
						    	A = Q / (U0 * deltaTm)
						    	DI = DE - 2 * G
						    	AT = 3.1415 * DE * L
						    	NT = A / AT
						    	VA = Np * FMA * (1 / NT) / (DA * 3.1415 * ((DI / 2)**2))
						    	ReA = DA * VA * DI / ViA
						    	if np.any(ReA < 2e4):
						    		F = 0.184 * (ReA)**(-1 / 5)
						    	if np.any(ReA >= 2e4):
						    		F = 0.184 * (ReA)**(-1 / 5)
						    	f = F / 8
						    	g = ReA - 1000
						    	h = (PrA**(2 / 3) - 1)
						    	ia = (F / 8)**(1 / 2)
						    	NuA = f * g * PrA / (1 + 12.7 * ia * h)
						    	HA = (4200 *(1.35+0.02*tm)*(VA**0.8))/((DI*(1000/1))**0.2) 
						    	Pt = 1.25 * DE
						    	if L == 1.83:
						    		NoBaffles = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29 ,30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40]
						    	if L == 2.44:
						    		NoBaffles = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29 ,30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40]
						    	if L == 3.66:
						    		NoBaffles = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29 ,30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40]
						    	if L == 4.88:
						    		NoBaffles = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29 ,30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40]
						    	if L == 6.10:
						    		NoBaffles = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29 ,30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40]
						    	if L == 7.32:
						    		NoBaffles = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29 , 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40]
						    		for NoBaffles in NoBaffles:
						    			lB = L / NoBaffles
						    			if Arreglo == 'Triangular':
						    				K1 = (0.0007 * ((Np)**4)) - (0.0117 * ((Np)**3)) + (0.0689 * ((Np)**2)) - (0.2054 * Np) + 0.4664
						    				n1=(0.0053 * (Np**2)) + (0.0285 * Np) + 2.1128
						    			if Arreglo == 'Cuadrado':
						    				K1 = (0.0017 * ((Np)**4)) - (0.0289 * ((Np)**3)) + (0.1633 * ((Np)**2)) - (0.372 * Np) + 0.4509
						    				n1 = (0.0734 * ((Np)**3)) - (0.0044 * ((Np)**4)) - (0.3922 * ((Np)**2)) + (0.183 * Np) + 1.7173
						    			DB = DE * (NT / K1)**(1 / n1)
						    			if ClaroCoraza == 'Pull-through floating head':
						    				Claro = 0.009375 * DB + 0.0856625
						    			if ClaroCoraza == 'Split-ring floating head':
						    				Claro = 0.02777 * DB + 0.0444
						    			if ClaroCoraza == 'Outside packed head':
						    				Claro = 0.0375
						    			if ClaroCoraza == 'Fixen and U-tube':
						    				Claro = 0.01 * DB + 0.008
						    			Ds = DB + Claro
						    			As = ((Pt - DE) * Ds * lB) / Pt
						    			VAD = FMAD / (As * DAD)
						    			if Arreglo == 'Triangular':
						    				de = (1.1 / DE) * ((Pt**2) - 0.917 * (DE**2))
						    			if Arreglo == 'Cuadrado':
						    				de = (1.27 / DE) * ((Pt**2) - 0.785 * (DE**2))
						    			ReAD = DAD * VAD * de / ViAD
						    			if CB == 0.15:
						    				jhCoraza = 0.4606 * ((ReAD)**-0.461)
						    			if CB == 0.25:
						    				jhCoraza = 0.4294 * ((ReAD)**-0.466)
						    			if CB == 0.35:
						    				jhCoraza = 0.3519 * ((ReAD)**-0.462)
						    			if CB == 0.45:
						    				jhCoraza = 0.2941 * ((ReAD)**-0.453)
						    			NuAD = (jhCoraza) * (ReAD) * (PrAD**(1 / 3))
						    			hAD = NuAD * KAD / de
						    			k = 60.5
						    			Ua = (1 / HA) * (DE / DI)
						    			Ub = 1 / hAD
						    			Uc = DE * (math.log(DE / DI)) / (2 * k)
						    			jfTubos = 0.0474 * (ReA**-0.248)
						    			deltaPTubos = Np * (8 * jfTubos * (L/DI) + 2.5) * (DA * (VA**2) / 2)
						    			if CB == 0.15:
						    				jfCoraza = 0.327 * ((ReAD)**-0.169) 
						    			if CB == 0.25:
						    				jfCoraza = 0.2263 * ((ReAD)**-0.163)
						    			if CB == 0.35:
						    				jfCoraza = 0.2187 * ((ReAD)**-0.177)
						    			if CB == 0.45:
						    				jfCoraza = 0.1801 * ((ReAD)**-0.179)
						    			deltaPCoraza = 8 * jfCoraza * (Ds / de) * (L / lB) * ((DAD * (VAD**2)) / 2)
						    			U1 = 1 / (Ua + Ub + Uc)
						    			Error = abs(U0 - U1) * 100 / U0
						    			if np.any(Error < 24):
						    				if np.any(deltaPTubos > 49900):
						    					if np.any(deltaPCoraza > 51000):
						    						if np.any(deltaPCoraza < 53080):
						    							if np.any(deltaPTubos < 51710):
						    								print('--------------------------------------------------------------------')
						    								print('Folio', i)
						    								print('U0', U0)
						    								print('El porcentaje de Error es del', Error, '%,')
						    								print('La caída de presión en la coraza es de', deltaPCoraza, 'Pa')
						    								print('La caída de presión en los tubos es de ', deltaPTubos, 'Pa')
						    								print('Diamétro exterior de', DE, 'm')
						    								print('Diámetro interior de', DI, 'm')
						    								print('Grosor de', G, 'm')
						    								print('Número de pasos', Np)
						    								print('Número de Baffles', NoBaffles-1)
						    								print('Longitud de', L)
						    								print('Corte de baffle', CB)
						    								print('Número de tubos', NT)
						    								print('Tipo de Claro en la coraza', ClaroCoraza)
						    								print('Arreglo de tubos', Arreglo)
						    								print('U1', U1)
						    								print('--------------------------------------------------------------------')
						    			i = i + 1
msvcrt.getch()