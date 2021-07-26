import kivy
kivy.require('1.0.7')

from kivy.app import App
from kivy.uix.gridlayout import GridLayout
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.label import Label
from kivy.uix.textinput import TextInput
from kivy.uix.button import Button
from kivy.config import Config

import configparser
import shutil
import os.path

configpath = "tmp.ini"

if not(os.path.isfile(configpath)):
    shutil.copy2('config.ini', 'tmp.ini')

Config.set('graphics', 'width', '200')
Config.set('graphics', 'height', '400')

config = configparser.ConfigParser()
config.read(configpath)

class IniForm(BoxLayout):

    def __init__(self, **kwargs):
        super(IniForm, self).__init__(**kwargs)
        self.orientation = 'vertical'

    def get_k(self):
        return config.get('Prop', 'k')

    def get_name(self):
        return config.get('Prop', 'name')

    def get_mMol(self):
        return config.get('Prop', 'mMol')
	
    def get_r(self):
	return config.get('Prop', 'r')

    def get_den_o(self):
	return config.get('Prop', 'den_o')

    def get_den_f(self):
	return config.get('Prop', 'den_f')

    def get_matl(self):
	return config.get('Nozzle', 'matl')

    def get_thk(self):
	return config.get('Nozzle', 'thk')

    def get_rho(self):
	return config.get('Nozzle', 'rho')

    def get_per(self):
        return config.get('Nozzle', 'per')

    def get_tp(self):
        return config.get('Nozzle', 'tp')
  
    def get_alpha(self):
	return config.get('Nozzle', 'alpha')

    def get_Tc(self):
        return config.get('Chamber', 'Tc')

    def get_Pc(self):
        return config.get('Chamber', 'Pc')

    def get_thetaconver(self):
	return config.get('Chamber', 'thetaconver')

    def write(self, name, k, mMol, r, den_o, den_f, per, matl, thk, tp, alpha, Tc, Pc, thetaconver):
        config.set('Prop', 'name', name.text)
        config.set('Prop', 'k', k.text)
        config.set('Prop', 'mMol', mMol.text)
	config.set('Prop', 'r', r.text)
	config.set('Prop', 'den_o', den_o.text)
	config.set('Prop', 'den_f', den_f.text)
        config.set('Nozzle', 'per', per.text)
        config.set('Nozzle', 'matl', matl.text)
	config.set('Nozzle', 'thk', thk.text)
	config.set('Nozzle', 'alpha', alpha.text)
        config.set('Nozzle', 'tp', tp.text)
        config.set('Chamber', 'Tc', Tc.text)
        config.set('Chamber', 'Pc', Pc.text)
	config.set('Chamber', 'tentaconver', thetaconver.text)

        with open(configpath, 'w') as configfile:
            config.write(configfile)
            

class NOZZLEDESIGNApp(App):
    def build(self):
        return IniForm()


if __name__ == '__main__':
    TDUApp().run()
