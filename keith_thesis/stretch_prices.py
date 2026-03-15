import pandas as pd
import os

file_in    = "year_lmp.csv"
output_dir = "stretched_lmp_cases"
col        = "LMP"

stretch_factors = [1.25, 1.5, 1.75, 2.0]

os.makedirs(output_dir, exist_ok=True)

for alpha in stretch_factors:
    df = pd.read_csv(file_in)
    prices = df[col].copy()

    p10 = prices.quantile(0.10)
    p90 = prices.quantile(0.90)

    new_prices = prices.copy()

    top_mask = prices >= p90
    new_prices.loc[top_mask] = p90 + alpha * (prices.loc[top_mask] - p90)

    bottom_mask = prices <= p10
    new_prices.loc[bottom_mask] = p10 - alpha * (p10 - prices.loc[bottom_mask])

    df_out      = df.copy()
    df_out[col] = new_prices

    yearly_file = os.path.join(output_dir, f"lmp_stretched_{alpha}.csv")
    df_out.to_csv(yearly_file, index=False)
    print(f"Saved {yearly_file}")


    df_avg = df_out.copy()
    df_avg["hour_of_day"] = [i % 24 for i in range(len(df_avg))]
    avg_day = df_avg.groupby("hour_of_day", as_index=False)["LMP"].mean()


    avg_file = os.path.join(output_dir, f"avg_day_lmp_{alpha}.csv")
    avg_day.to_csv(avg_file, index=False)
    print(f"Saved {avg_file}")