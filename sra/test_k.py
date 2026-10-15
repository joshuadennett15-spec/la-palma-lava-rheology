import pandas as pd
import numpy as np

df = pd.read_csv('sheet2_out.csv')
g=9.81
rho=2407.5
vent_elevs = df.groupby('Eruption')['Elevation (m)'].max().to_dict()
df['E_vent'] = df['Eruption'].map(vent_elevs)

for k in [0.001, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03]:
    df['Q_local'] = df['Q (m³/s)'] * np.exp(-k * (df['E_vent'] - df['Elevation (m)']))
    df['n'] = 3*(1+df['Height_Mean']/df['Width_Mean'])**2
    df['eta'] = (rho*g*df['Height_Mean']**3*df['Width_Mean']*np.sin(np.radians(df['Slope_Mean'])))/(df['n']*df['Q_local'])
    
    # print the first and last viscosity for both
    sanjuan = df[df['Eruption']=='1949_San_Juan']
    teneguia = df[df['Eruption']=='1971_Teneguia']
    
    sj_first = sanjuan.iloc[0]['eta']
    sj_last = sanjuan.iloc[-1]['eta']
    t_first = teneguia.iloc[0]['eta']
    t_last = teneguia.iloc[-1]['eta']
    print(f"k={k}:")
    print(f"  San Juan start: {sj_first:.2e}, end: {sj_last:.2e} (ratio {sj_last/sj_first:.2f})")
    print(f"  Teneguia start: {t_first:.2e}, end: {t_last:.2e} (ratio {t_last/t_first:.2f})")
