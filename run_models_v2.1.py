import numpy as np
import subprocess

import time

time_start=time.time()


para_list = []
# for M in [0.12,0.15,0.20,0.28]:
#     for ϵ in [0.0, 0.5e-13,1e-13,2e-13]:#
#         para_list.append({"M":M, "ϵ":ϵ,"m":1000.0})

eps_m_curve=np.log10(np.load("eps-m-curve.npy"))

from scipy.interpolate import interp1d
def estimate_ϵ(M,m):
    iM=round(np.interp(M,[0.12,0.15,0.20,0.28],[0,1,2,3]))
    lgm=np.log10(m)

    f = interp1d(eps_m_curve[iM][0,:], eps_m_curve[iM][3,:],
    kind='linear',
    bounds_error=False,
    fill_value='extrapolate')
    return 10**f(lgm)


for M in [0.12,0.15,0.20,0.28]:
    # for m in np.geomspace(1,1e4,40):
    for m in np.geomspace(150,650,10):
        ϵ_est=estimate_ϵ(M,m)
        for ϵ in [ϵ_est]:#
            para_list.append({"M":M, "ϵ":ϵ,"m":m})

# for M in [0.12,0.15,0.20,0.28]:
#     for m in [800.0]: #np.geomspace(1,1e4,20):
#         ϵ_est=estimate_ϵ(M,m)
#         for ϵ in [0.0, 0.7*ϵ_est, ϵ_est, 1.4*ϵ_est ]:#
#             para_list.append({"M":M, "ϵ":ϵ,"m":m})

# for M in [0.12,0.15,0.20,0.28,0.32]:
#     for ϵ in [0.0]:#
#         for m in [1.0]:
#             para_list.append({"M":M, "ϵ":ϵ,"m":m})

print("length of para_list=",len(para_list))


# store results
results = []

# read template inlist
with open("inlist_template_v2", "r") as f:
    template = f.read()

for para in para_list:

    M=para["M"]
    ϵ=para["ϵ"]
    m=para["m"]
    print(f"Running:", para)

    # create inlist_project
    inlist_text = template.replace("MASS_MSUN", str(M))
    inlist_text = inlist_text.replace("DP_mass_in_eV", str(m))
    inlist_text = inlist_text.replace("DP_epsilon", str(ϵ))
    

    with open("inlist_project", "w") as f:
        f.write(inlist_text)

    # run MESA
    with open(f"run.log", "w") as logfile:
        subprocess.run(
            ["./rn"],
            stdout=logfile,
            stderr=logfile
        )

    #save final profile
    data = np.genfromtxt(
        "final_profile.data",
        comments="#",
        skip_header=5,
        names=True
    )
    np.save(f"result/M{M:.2f}-m{m:.2f}-e{ϵ/1e-15:.1f}-profile.npy",data)

    # read history.data
    data = np.genfromtxt(
        "LOGS/history.data",
        comments="#",
        skip_header=5,
        names=True
    )

    np.save(f"result/M{M:.2f}-m{m:.2f}-e{ϵ/1e-15:.1f}.npy",data)


    

    # extract final radius
    final_R = 10**(data['log_R'][-1])

    # save result
    results.append([m, ϵ, M, final_R])



# convert to numpy array
results = np.array(results)

# save output
np.savetxt(
    "output.dat",
    results,
    header="DP_m[eV] DP_epsilon  M[Msol]    R[Rsol]"
)



print("Finished. Results saved in output.dat")
Δt=time.time()-time_start
print("time used:",Δt)