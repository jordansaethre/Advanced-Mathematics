import numpy as np
import pandas as pd
import plotly.graph_objects as go

# parameters
deltaxy1 = 1/100 # delta x = delta y
deltaxy2 = 1/200
deltaxy3 = 1/400
step1 = 200
step2 = 400
step3 = 800
deltat = 0.012
iterates1 = 1
iterates2 = 10
s1 = (deltat/(deltaxy1**2)) # mu_x = mu_y
s2 = (deltat/(deltaxy2**2))
s3 = (deltat/(deltaxy3**2))

# initial condition

x = lambda theta: (0.5+ 0.1*np.sin(6*theta))*np.cos(theta)

y = lambda theta: (0.5+ 0.1*np.sin(6*theta))*np.sin(theta)


def initial_mat(x,y,n):
    initial=np.zeros(shape=(n,n))
    inter=np.linspace(-1,1,n)
    for i in range(0,len(inter)):
        for j in range(0,len(inter)):
            ri=np.sqrt(inter[i]**2 +inter[j]**2)
            r=np.sqrt(x(np.arccos(inter[i]/ri))**2 + y(np.arccos(inter[i]/ri))**2)
            if abs(ri)>abs(r):
                initial[i,j]=0
            else:
                initial[i,j]=1
    return(initial)

# tridiagonal matrix
def tridag(L,M,U,k1=-1,k2=0,k3=1):
    return(np.diag(L,k1)+np.diag(M,k2)+np.diag(U,k3))

# one adi xy sweep
def adi_method(s, deltaxy, step, initial):
	# diagonals
	lower=[-s/2 for i in range(0,len(initial)-1)]
	main=[1+s for i in range(0,len(initial))]
	upper=[-s/2 for i in range(0,len(initial)-1)]

	# tridiagonal matrix
	M = tridag(lower,main,upper)
	# calculate inverse
	Minv = np.linalg.inv(M)

	initial = np.transpose(initial)
	# x sweep
	for i in range(0, step):
		initial[i] = np.matmul(Minv,np.transpose(initial[i]))
		i = i + 1
	# y sweep	
	initial = np.transpose(initial)
	for j in range(0, step):
		initial[j] = np.matmul(Minv,np.transpose(initial[j]))
		j = j + 1

	initial = np.transpose(initial)

	return initial

def adi_loop(s, deltaxy, step, initial, iterates):
	solution = adi_method(s, deltaxy, step, initial)
	for n in range(0,iterates):
		adi_method(s, deltaxy, step, solution)
		n = n + 1
	return solution

def plot_surface(solution):
	surface_df = pd.DataFrame(solution)
	surface_df.to_csv('surface_data.csv', index = False)
	surface_data = pd.read_csv('surface_data.csv')
	fig = go.Figure(data=[go.Surface(z=surface_data.values, colorscale='Viridis')])
	fig.update_traces(contours_z=dict(show=True, usecolormap=True, highlightcolor="limegreen", project_z=True))
	fig.update_layout(title='Solution', autosize=False, width=500, height=500,margin=dict(l=65, r=50, b=65, t=90), xaxis = dict(visible = False))
	fig.show()


# deltax = deltay = 1/100
# # initial
# plot_surface(initial_mat(x,y,step1))
# ## after 1 iterations
# plot_surface(adi_loop(s1, deltaxy1, step1, initial_mat(x,y,step1), iterates1))
# # after 10 iterations
# plot_surface(adi_loop(s1, deltaxy1, step1, initial_mat(x,y,step1), iterates2))

# deltax = deltay = 1/200
# # initial
# plot_surface(initial_mat(x,y,step2))
## after 1 iterations
# plot_surface(adi_loop(s2, deltaxy2, step2, initial_mat(x,y,step2), iterates1))
# after 10 iterations
# plot_surface(adi_loop(s2, deltaxy2, step2, initial_mat(x,y,step2), iterates2))

# # deltax = deltay = 1/400
# ## initial
# plot_surface(initial_mat(x,y,step3))
# ## after 1 iterations
# plot_surface(adi_loop(s3, deltaxy3, step3, initial_mat(x,y,step3), iterates1))
# ## after 10 iterations
plot_surface(adi_loop(s3, deltaxy3, step3, initial_mat(x,y,step3), iterates2))