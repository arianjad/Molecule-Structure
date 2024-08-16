import multiprocessing as mp
import numpy as np
from qutip import *

def prep_sim(matrix_info,state_groups,qops,qd,s_X=1,s_Z=1,one_photon=0,e_decohere_rate = 10*2*np.pi, g_decohere_rate=0.05*2*np.pi, decay_rate = 1/0.02,vbr_scaling=None,time = np.linspace(0,4,5000),big_decohere=False,pol='XZ'):

    H0,Hoff_Z,Hoff_X,BR0,BRp,BRm,rho0 = matrix_info
    e_states = state_groups[0]
    g_states = [g for group in state_groups[1:4] for g in group]

    Rabi_Z = decay_rate*np.sqrt(1/2*s_Z)/2
    Rabi_X =  decay_rate*np.sqrt(1/2*s_X)/2

    if vbr_scaling is not None:
        Rabi_Z*=np.sqrt(vbr_scaling)
        Rabi_X*=np.sqrt(vbr_scaling)

    g_decohere = [np.sqrt(g_decohere_rate)*qops[i][i] for i in range(len(g_states))]
    e_decohere = [np.sqrt(e_decohere_rate)*qops[i][i] for i in range(len(e_states))]
    decay = [np.sqrt(decay_rate)*(Qobj(BR.T)) for BR in [BR0,BRp,BRm]]
    if big_decohere == False:
        collapse=[*g_decohere,*e_decohere,*decay]
    else:
        collapse=[*g_decohere,*e_decohere,*decay,*[[qops[i][i]*np.sqrt(10),on_off2] for i in range(len(qops))]]
    size = len(qops)
    Hd = np.zeros((size,size))
    for i in range(len(state_groups[0])):
        Hd[i,i]-= one_photon
    if pol=='XZ' or pol=='ZX':
        H=[Qobj(H0+Hd)*2*np.pi,[Rabi_Z*Qobj(Hoff_Z + Hoff_Z.T),on_off],[Rabi_X*Qobj(Hoff_X),exp_onoff],[Rabi_X*Qobj(Hoff_X.T),expconj_onoff]]
    else:
        H=[Qobj(H0+Hd)*2*np.pi,[Rabi_Z*Qobj(Hoff_Z + Hoff_Z.T),on_off],[Rabi_X*Qobj(Hoff_Z),exp_onoff],[Rabi_X*Qobj(Hoff_Z.T),expconj_onoff]]
    return H, rho0, time, collapse, state_groups

def mp_cpt_sim(two_photon_array,worker_func,H,rho0,time,collapse,state_groups,tf=3):
    pool = mp.Pool(mp.cpu_count())
    results = pool.starmap(worker_func, [(two_photon,H,rho0,time,collapse,state_groups,tf) for two_photon in two_photon_array])
    pool.close()
    pool.join()
    result_array = np.array(results)
    return result_array

def prep_sim_1ph(one_photon_array,matrix_info,state_groups,qops,qd,s_X=1,s_Z=1,e_decohere_rate = 10*2*np.pi, g_decohere_rate=0.05*2*np.pi, decay_rate = 1/0.02,vbr_scaling=None,time = np.linspace(0,4,5000),big_decohere=False,pol='XZ'):

    H0,Hoff_Z,Hoff_X,BR0,BRp,BRm,rho0 = matrix_info
    e_states = state_groups[0]
    g_states = [g for group in state_groups[1:4] for g in group]

    Rabi_Z = decay_rate*np.sqrt(1/2*s_Z)/2
    Rabi_X =  decay_rate*np.sqrt(1/2*s_X)/2

    if vbr_scaling is not None:
        Rabi_Z*=np.sqrt(vbr_scaling)
        Rabi_X*=np.sqrt(vbr_scaling)

    g_decohere = [np.sqrt(g_decohere_rate)*qops[i][i] for i in range(len(g_states))]
    e_decohere = [np.sqrt(e_decohere_rate)*qops[i][i] for i in range(len(e_states))]
    decay = [np.sqrt(decay_rate)*(Qobj(BR.T)) for BR in [BR0,BRp,BRm]]
    if big_decohere == False:
        collapse=[*g_decohere,*e_decohere,*decay]
    else:
        collapse=[*g_decohere,*e_decohere,*decay,*[[qops[i][i]*np.sqrt(10),on_off2] for i in range(len(qops))]]

    size = len(qops)
    H_all = []
    for one_photon in one_photon_array:
        Hd = np.zeros((size,size))
        for i in range(len(state_groups[0])):
            Hd[i,i]-= one_photon
        if pol=='XZ' or pol=='ZX':
            H=[Qobj(H0+Hd)*2*np.pi,[Rabi_Z*Qobj(Hoff_Z + Hoff_Z.T),on_off],[Rabi_X*Qobj(Hoff_X),exp_onoff],[Rabi_X*Qobj(Hoff_X.T),expconj_onoff]]
        else:
            H=[Qobj(H0+Hd)*2*np.pi,[Rabi_Z*Qobj(Hoff_Z + Hoff_Z.T),on_off],[Rabi_X*Qobj(Hoff_Z),exp_onoff],[Rabi_X*Qobj(Hoff_Z.T),expconj_onoff]]
        H_all.append(H)
    return H_all, rho0, time, collapse, state_groups

def mp_cpt_sim_1ph(one_photon_array,two_photon,worker_func,H,rho0,time,collapse,state_groups,tf=3):
    pool = mp.Pool(mp.cpu_count())
    results = pool.starmap(worker_func, [(two_photon,H_,rho0,time,collapse,state_groups,tf) for H_ in H])
    pool.close()
    pool.join()
    result_array = np.array(results)
    return result_array

def cpt_sim(two_photon, H,rho0,time,collapse,state_groups,tf):
    two_photon = two_photon*2*np.pi
    result = mesolve(H,rho0,time,collapse,[sum(states) for states in state_groups],args={'on':0,'off':tf,'on2':3.2,'off2':4,'w':two_photon})
    return result.expect

def prep_sim_1d_weighted(matrix_info,state_groups,qops,qd,s_X=1,s_Z=1,one_photon=0,e_decohere_rate = 10*2*np.pi, g_decohere_rate=0.05*2*np.pi, decay_rate = 1/0.02,vbr_scaling=None,time = np.linspace(0,4,5000),big_decohere=False):

    H0,Hoff_Z,Hoff_X,BR0,BRp,BRm,rho0 = matrix_info
    e_states = state_groups[0]
    g_states = [g for group in state_groups[1:4] for g in group]

    Rabi_Z = decay_rate*np.sqrt(1/2*s_Z)/2
    Rabi_X =  decay_rate*np.sqrt(1/2*s_X)/2

    if vbr_scaling is not None:
        Rabi_Z*=np.sqrt(vbr_scaling)
        Rabi_X*=np.sqrt(vbr_scaling)

    g_decohere = [np.sqrt(g_decohere_rate)*qops[i][i] for i in range(len(g_states))]
    e_decohere = [np.sqrt(e_decohere_rate)*qops[i][i] for i in range(len(e_states))]
    decay = [np.sqrt(decay_rate)*(Qobj(BR.T)) for BR in [BR0,BRp,BRm]]
    if big_decohere == False:
        collapse=[*g_decohere,*e_decohere,*decay]
    else:
        collapse=[*g_decohere,*e_decohere,*decay,*[[qops[i][i]*np.sqrt(10),on_off2] for i in range(len(qops))]]

    size = len(qops)
    Hd = np.zeros((size,size))
    for i in range(len(state_groups[0])):
        Hd[i,i]-= one_photon

    H=Qobj(H0+Hd)*2*np.pi
    Hall = []
    Hinfo = [H,Hoff_Z,Hoff_X,Rabi_X,Rabi_Z]
    # for s_ in s_array:
    #     Hs = [H,[decay_rate*np.sqrt(1/2*s_)/2*Qobj(Hoff_Z + Hoff_Z.T),on_off],[decay_rate*np.sqrt(1/2*s_)/2*Qobj(Hoff_X),exp_onoff],[decay_rate*np.sqrt(1/2*s_)/2*Qobj(Hoff_X.T),expconj_onoff]]
    #     Hall.append(Hs)

    return Hinfo, rho0, time, collapse, state_groups

def mp_cpt_sim_1d_weighted(two_photon_array,s_weight_array,worker_func,Hinfo,rho0,time,collapse,state_groups,tf=3,pol='XZ'):
    pool = mp.Pool(mp.cpu_count())
    results = pool.starmap(worker_func, [(two_photon,s_weight,Hinfo,rho0,time,collapse,state_groups,tf,pol) for two_photon,s_weight in zip(two_photon_array,s_weight_array)])
    pool.close()
    pool.join()
    result_array = np.array(results)
    return result_array


def cpt_sim_1d_weighted(two_photon, s_weight ,Hinfo,rho0,time,collapse,state_groups,tf,pol):
    two_photon = two_photon*2*np.pi

    H0, Hoff_Z, Hoff_X,Rabi_X,Rabi_Z = Hinfo
    if pol=='XZ':
        H = [H0,[Rabi_Z*Qobj(Hoff_Z + Hoff_Z.T),on_off],[Rabi_X*np.sqrt(s_weight)*Qobj(Hoff_X),exp_onoff],[Rabi_X*np.sqrt(s_weight)*Qobj(Hoff_X.T),expconj_onoff]]
    else:
        H = [H0,[Rabi_Z*Qobj(Hoff_Z + Hoff_Z.T),on_off],[Rabi_X*np.sqrt(s_weight)*Qobj(Hoff_Z),exp_onoff],[Rabi_X*np.sqrt(s_weight)*Qobj(Hoff_Z.T),expconj_onoff]]
    result = mesolve(H,rho0,time,collapse,[sum(states) for states in state_groups],args={'on':0,'off':tf,'w':two_photon})
    return result.expect



def prep_sim2d(s_array,matrix_info,state_groups,qops,qd,one_photon=0,e_decohere_rate = 10*2*np.pi, g_decohere_rate=0.05*2*np.pi, decay_rate = 1/0.02,vbr_scaling=None,time = np.linspace(0,4,5000),big_decohere=False):

    H0,Hoff_Z,Hoff_X,BR0,BRp,BRm,rho0 = matrix_info
    e_states = state_groups[0]
    g_states = [g for group in state_groups[1:4] for g in group]

    Rabi_Z = decay_rate*np.sqrt(1/2*s_Z)/2
    Rabi_X =  decay_rate*np.sqrt(1/2*s_X)/2

    if vbr_scaling is not None:
        Rabi_Z*=np.sqrt(vbr_scaling)
        Rabi_X*=np.sqrt(vbr_scaling)

    g_decohere = [np.sqrt(g_decohere_rate)*qops[i][i] for i in range(len(g_states))]
    e_decohere = [np.sqrt(e_decohere_rate)*qops[i][i] for i in range(len(e_states))]
    decay = [np.sqrt(decay_rate)*(Qobj(BR.T)) for BR in [BR0,BRp,BRm]]
    if big_decohere == False:
        collapse=[*g_decohere,*e_decohere,*decay]
    else:
        collapse=[*g_decohere,*e_decohere,*decay,*[[qops[i][i]*np.sqrt(10),on_off2] for i in range(len(qops))]]

    size = len(qops)
    Hd = np.zeros((size,size))
    for i in range(len(state_groups[0])):
        Hd[i,i]-= one_photon

    H=Qobj(H0+Hd)*2*np.pi
    Hall = []
    Hinfo = [H,Hoff_Z,Hoff_X,decay_rate]
    # for s_ in s_array:
    #     Hs = [H,[decay_rate*np.sqrt(1/2*s_)/2*Qobj(Hoff_Z + Hoff_Z.T),on_off],[decay_rate*np.sqrt(1/2*s_)/2*Qobj(Hoff_X),exp_onoff],[decay_rate*np.sqrt(1/2*s_)/2*Qobj(Hoff_X.T),expconj_onoff]]
    #     Hall.append(Hs)

    return Hinfo, rho0, time, collapse, state_groups

def mp_cpt_sim_2d(two_photon_array,s_array,worker_func,Hinfo,rho0,time,collapse,state_groups,tf=3,pol='XZ'):
    mesh = np.array(np.meshgrid(two_photon_array,s_array,indexing='ij'))
    pairs = np.stack((mesh[0].ravel(),mesh[1].ravel()),axis=1)
    pool = mp.Pool(mp.cpu_count())
    results = pool.starmap(worker_func, [(two_photon,s_,Hinfo,rho0,time,collapse,state_groups,tf,pol) for two_photon,s_ in pairs])
    pool.close()
    pool.join()
    result_array = np.array(results)
    return result_array


# def cpt_sim_2d(two_photon, s ,Hinfo,rho0,time,collapse,state_groups,tf,pol):
#     two_photon = two_photon*2*np.pi
#     H0, Hoff_Z, Hoff_X, decay_rate = Hinfo
#     if pol=='XZ':
#         H = [H0,[decay_rate*np.sqrt(1/2*s)/2*Qobj(Hoff_Z + Hoff_Z.T),on_off],[decay_rate*np.sqrt(1/2*s)/2*Qobj(Hoff_X),exp_onoff],[decay_rate*np.sqrt(1/2*s)/2*Qobj(Hoff_X.T),expconj_onoff]]
#     else:
#         H = [H0,[decay_rate*np.sqrt(1/2*s)/2*Qobj(Hoff_Z + Hoff_Z.T),on_off],[decay_rate*np.sqrt(1/2*s)/2*Qobj(Hoff_Z),exp_onoff],[decay_rate*np.sqrt(1/2*s)/2*Qobj(Hoff_Z.T),expconj_onoff]]
#     result = mesolve(H,rho0,time,collapse,[sum(states) for states in state_groups],args={'on':0,'off':tf,'on2':3.2,'off2':4,'w':two_photon})
#     return result.expect


def cpt_sim_2d(two_photon, s ,Hinfo,rho0,time,collapse,state_groups,tf,pol):
    two_photon = two_photon*2*np.pi
    H0, Hoff_Z, Hoff_X, decay_rate = Hinfo
    if pol=='XZ':
        H = [H0,[decay_rate*np.sqrt(1/2*s)/2*Qobj(Hoff_Z + Hoff_Z.T),on_off],[decay_rate*np.sqrt(1/2*s)/2*Qobj(Hoff_X),exp_onoff],[decay_rate*np.sqrt(1/2*s)/2*Qobj(Hoff_X.T),expconj_onoff]]
    else:
        H = [H0,[decay_rate*np.sqrt(1/2*s)/2*Qobj(Hoff_Z + Hoff_Z.T),on_off],[decay_rate*np.sqrt(1/2*s)/2*Qobj(Hoff_Z),exp_onoff],[decay_rate*np.sqrt(1/2*s)/2*Qobj(Hoff_Z.T),expconj_onoff]]
    result = mesolve(H,rho0,time,collapse,[sum(states) for states in state_groups],args={'on':0,'off':tf,'on2':3.2,'off2':4,'w':two_photon})
    return result.expect

# generate n level basis
def q_basis(n):
    q_basis = []
    for i in range(n):
        q_basis.append(basis(n,i))
    return q_basis

def q_dict(q_basis, q_str):
    q_dict = {qs: qb for qs,qb in zip(q_str, q_basis)}
    return q_dict

def q_ops(q_basis):
    n = len(q_basis)
    q_ops = np.zeros((n,n)).tolist()
    for i in range(n):
        for j in range(n):
            q_ops[i][j] = q_basis[i]*q_basis[j].dag()
    return q_ops

def exp_osc(t,args):
    w = args['w']
    return np.exp(1j*w*t)

def exp_conj(t,args):
    w = args['w']
    return np.exp(-1j*w*t)

def linear_ramp(t,args):
    m = args['slope']
    b = args['start']
    stop = args['stop']
    y = m*t+b
    if y>stop:
        return 0
    else:
        return y

def on_off(t,args):
    on = args['on']
    off = args['off']
    if t>on and t<off:
        return 1
    else:
        return 0

def on_off2(t,args):
    on = args['on2']
    off = args['off2']
    if t>on and t<off:
        return 1
    else:
        return 0

def exp_onoff(t,args):
    w = args['w']
    on = args['on']
    off = args['off']
    exp = np.exp(1j*w*t)
    if t>on and t<off:
        return 1*exp
    else:
        return 0

def expconj_onoff(t,args):
    w = args['w']
    on = args['on']
    off = args['off']
    exp = np.exp(-1j*w*t)
    if t>on and t<off:
        return 1*exp
    else:
        return 0
