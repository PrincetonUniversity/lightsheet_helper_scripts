#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 18:05:48 2019

@author: wanglab
"""

import os, time, smtplib, ssl
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

#this is run on the ls acquisition computer

def send_email(img, receiver_email):
    file = os.stat(img)
    mod = time.ctime(file.st_mtime)
    
    if receiver_email == "": receiver_email = "zmd@princeton.edu"
    
    print("\nlast modified time: ", mod)
    
    if time.time() - file.st_mtime > 180: #if its not been modified for over 3 minutes 
        
        smtp_server = "smtp.gmail.com"
        port = 587  # For starttls
        sender_email = "pni.lvbt.lightsheet@gmail.com"
        password = "stephanisawesome"
        msg = MIMEMultipart()
        msg["Subject"] = "Light-sheet acquisition complete!"
        body = "Last modified time of file: {}".format(mod)
        msg.attach(MIMEText(body, 'plain'))
        context = ssl.create_default_context()
        with smtplib.SMTP(smtp_server, port) as server:
            server.ehlo()  # Can be omitted
            server.starttls(context=context)
            server.ehlo()  # Can be omitted
            server.login(sender_email, password)
            text = msg.as_string()
            server.sendmail(sender_email, receiver_email, text)
        print ("\nsent email :D \n")
    
    else:
        print("acquisition not complete! \n")

if __name__ == "__main__":
    
    img = "G:\\190307_20170130_tp_bl6_sim_rpv_01_4x_647_008na_1hfds_z7d5um_150msec_10povlp_14-58-10"
        
    receiver_email = "zmd@princeton.edu"
    
    send_email(img, receiver_email)