sudo ln -fs /etc/nginx/sites-available/indra_network_search /etc/nginx/sites-enabled/
sudo rm /etc/nginx/sites-enabled/default
sudo systemctl start nginx

