import time
from datetime import datetime, timedelta
from concurrent.futures import ThreadPoolExecutor, as_completed
from run_forecast import create_forecast_job

# Define the time range
start_date = datetime(2012, 12, 7)
end_date = datetime(2024, 11, 1)
interval = timedelta(days=1)
max_concurrent_jobs = 9

def generate_date_ranges(start, end, interval):
    current = start
    while current < end:
        yield current, min(current + interval, end)
        current += interval

if __name__ == "__main__":
    # Using ThreadPoolExecutor to manage up to 5 concurrent jobs
    with ThreadPoolExecutor(max_workers=max_concurrent_jobs) as executor:
        futures = []
        
        # Generate jobs for each 10-day interval
        for start, end in generate_date_ranges(start_date, end_date, interval):
            futures.append(executor.submit(create_forecast_job, start, end))
        
        # Wait for all job1s to complete
        for future in as_completed(futures):
            unique_id, output_json, log_file = future.result()
            print(f"Job {unique_id} completed. Job JSON: {output_json}, Log: {log_file}")

