import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, cross_val_predict, KFold
from sklearn.feature_selection import VarianceThreshold, RFE
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.preprocessing import StandardScaler
import joblib

def build_qsar_model():
    print("Memuat data deskriptor...")
    df = pd.read_csv('h1_antagonists_descriptors.csv')
    
    # Hapus baris yang mengandung NaN
    df = df.dropna()
    print(f"Data setelah hapus NaN: {len(df)} baris")
    
    # Pisahkan X (deskriptor) dan y (aktivitas)
    X = df.drop(columns=['canonical_smiles', 'pIC50'])
    y = df['pIC50']
    
    # 1. Hapus deskriptor dengan variansi rendah (mendekati konstan)
    print("Menghapus fitur dengan variansi rendah...")
    vt = VarianceThreshold(threshold=0.01)
    X_vt = vt.fit_transform(X)
    X = pd.DataFrame(X_vt, columns=X.columns[vt.get_support()])
    print(f"Jumlah fitur tersisa: {X.shape[1]}")
    
    # 2. Hapus fitur yang sangat berkorelasi (Multikolinearitas > 0.85)
    # Ini sangat penting untuk Multiple Linear Regression (MLR)
    print("Menghapus fitur dengan korelasi tinggi...")
    corr_matrix = X.corr().abs()
    upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
    to_drop = [column for column in upper.columns if any(upper[column] > 0.85)]
    X = X.drop(columns=to_drop)
    print(f"Jumlah fitur tersisa setelah hapus korelasi: {X.shape[1]}")
    
    # 3. Pembagian Data (Train/Test Split 80:20)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    
    # Standarisasi (Skala Z-score)
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    X_train_scaled = pd.DataFrame(X_train_scaled, columns=X_train.columns)
    X_test_scaled = pd.DataFrame(X_test_scaled, columns=X_test.columns)
    
    # 4. Feature Selection dengan RFE (Recursive Feature Elimination)
    # Memilih 5 deskriptor terbaik untuk menghindari overfitting pada MLR
    print("Memilih 5 deskriptor terbaik dengan RFE untuk MLR...")
    estimator = LinearRegression()
    selector = RFE(estimator, n_features_to_select=5, step=1)
    selector = selector.fit(X_train_scaled, y_train)
    
    selected_features = X_train.columns[selector.support_].tolist()
    print(f"Fitur terpilih: {selected_features}")
    
    # 5. Membangun Model MLR
    X_train_final = X_train_scaled[selected_features]
    X_test_final = X_test_scaled[selected_features]
    
    mlr_model = LinearRegression()
    mlr_model.fit(X_train_final, y_train)
    
    # --- Validasi Internal ---
    # R2 Training
    y_train_pred = mlr_model.predict(X_train_final)
    r2_train = r2_score(y_train, y_train_pred)
    rmse_train = np.sqrt(mean_squared_error(y_train, y_train_pred))
    
    # Q2 (Cross-Validation R2) menggunakan 5-Fold
    cv = KFold(n_splits=5, shuffle=True, random_state=42)
    y_cv_pred = cross_val_predict(mlr_model, X_train_final, y_train, cv=cv)
    q2_cv = r2_score(y_train, y_cv_pred)
    
    # --- Validasi Eksternal ---
    y_test_pred = mlr_model.predict(X_test_final)
    r2_test = r2_score(y_test, y_test_pred)
    rmse_test = np.sqrt(mean_squared_error(y_test, y_test_pred))
    
    print("\n==================================")
    print("--- HASIL VALIDASI QSAR (MLR) ---")
    print("==================================")
    print(f"R2 Training    : {r2_train:.3f}")
    print(f"RMSE Training  : {rmse_train:.3f}")
    print(f"Q2 (5-Fold CV) : {q2_cv:.3f}")
    print(f"R2 Test (Ext)  : {r2_test:.3f}")
    print(f"RMSE Test      : {rmse_test:.3f}")
    print("==================================")
    
    # Tampilkan Persamaan MLR
    coef = mlr_model.coef_
    intercept = mlr_model.intercept_
    eq = f"pIC50 = {intercept:.3f}"
    for f, c in zip(selected_features, coef):
        sign = "+" if c >= 0 else "-"
        eq += f" {sign} {abs(c):.3f} * [{f}]"
    print("\nPersamaan QSAR MLR (Deskriptor terstandarisasi):")
    print(eq)
    
    # Menyimpan model
    joblib.dump(mlr_model, 'qsar_mlr_model.pkl')
    joblib.dump(scaler, 'qsar_scaler.pkl')
    joblib.dump(selected_features, 'qsar_selected_features.pkl')
    print("\nModel berhasil disimpan.")

if __name__ == "__main__":
    build_qsar_model()
