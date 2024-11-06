import json
import subprocess
from datetime import datetime, timedelta
import os

def create_forecast_job(start_date, end_date):
    # Load the JSON template
    with open('jobs/gfsa_template.json', 'r') as file:
        data = json.load(file)
    
    # Format start and end dates
    start_str = start_date.strftime("%Y-%m-%d_%H-%M-%S")
    end_str = end_date.strftime("%Y-%m-%d_%H-%M-%S")
    
    # Update JSON content
    data["start_utc"] = start_str
    data["end_utc"] = end_str
    
    # Generate filename based on start and end dates
    unique_name = f"{start_str}_to_{end_str}"
    output_json = f'jobs/{unique_name}.json'
    log_file = f'logs/{unique_name}.log'
    
    # Save modified JSON file
    with open(output_json, 'w') as file:
        json.dump(data, file, indent=4)
    
    # Run the forecast script
    with open(log_file, 'w') as logfile:
        subprocess.run(['./forecast.sh', output_json], stdout=logfile, stderr=logfile)

    return unique_name, output_json, log_file

if __name__ == "__main__":
    # Example usage for a 10-day period (modify as needed)
    start_date = datetime(2011, 1, 1)
    end_date = start_date + timedelta(days=10)
    create_forecast_job(start_date, end_date)

