import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import joblib

def plot_applicability_domain():
    print("Membuat plot Applicability Domain (Williams Plot)...")
    # Load model and data
    model = joblib.load('qsar_mlr_model.pkl')
    scaler = joblib.load('qsar_scaler.pkl')
    features = joblib.load('qsar_selected_features.pkl')
    
    df = pd.read_csv('h1_antagonists_descriptors.csv').dropna()
    X = df[features]
    y = df['pIC50']
    
    from sklearn.preprocessing import StandardScaler
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    y_pred = model.predict(X_scaled)
    residuals = y - y_pred
    
    # Standarisasi residual
    std_residuals = (residuals - np.mean(residuals)) / np.std(residuals)
    
    # Menghitung Leverage (Hat matrix)
    # H = X(X'X)^-1 X'
    X_mat = np.matrix(X_scaled)
    # Tambahkan intercept term (bias)
    X_mat = np.hstack([np.ones((X_mat.shape[0], 1)), X_mat])
    
    H = X_mat * np.linalg.pinv(X_mat.T * X_mat) * X_mat.T
    leverage = np.array(np.diag(H)).flatten()
    
    # Batas peringatan leverage (h*)
    k = len(features) # jumlah prediktor
    n = X_mat.shape[0] # jumlah sampel
    h_star = 3 * (k + 1) / n
    
    plt.figure(figsize=(10, 6))
    plt.scatter(leverage, std_residuals, alpha=0.7, edgecolor='k')
    plt.axhline(y=3, color='r', linestyle='--', label='±3 Standard Deviation')
    plt.axhline(y=-3, color='r', linestyle='--')
    plt.axvline(x=h_star, color='g', linestyle='--', label=f'Warning Leverage (h* = {h_star:.3f})')
    
    plt.xlabel('Leverage (h)', fontsize=12)
    plt.ylabel('Standardized Residuals', fontsize=12)
    plt.title('Williams Plot (Applicability Domain)', fontsize=14)
    plt.legend()
    
    plt.savefig('williams_plot.png', dpi=300, bbox_inches='tight')
    print(f"Plot Williams berhasil disimpan ke williams_plot.png")
    print(f"Batas Leverage (h*) = {h_star:.3f}")

if __name__ == "__main__":
    plot_applicability_domain()
