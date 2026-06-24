import numpy as np
import subprocess


para_list = []
M=1.0
# for m in np.geomspace(1,1e3,30):
# for m in np.geomspace(1e3,1e4,10):
for m in [0.2,0.4,0.8]:
    for ϵ in [1e-15]:#
        para_list.append({"M":M, "ϵ":ϵ,"m":m})
        
print("length of para_list=",len(para_list))


# store results
results = []

# read template inlist
with open("inlist_template_Sun", "r") as f:
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

    # read history.data
    data = np.genfromtxt(
        "LOGS/history.data",
        comments="#",
        skip_header=5,
        names=True
    )

    np.save(f"result/Sun-m{m}-e{ϵ/1e-15:.1f}.npy",data)

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