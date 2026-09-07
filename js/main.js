import { CHEM } from "./chem-content.js";
import { mountApp } from "./engine.js";
import { mountExtras } from "./chem-extras.js";
import { mountAccount } from "./account.js";
import { mountStats } from "./stats.js";
import { mountHome } from "./home.js";
import { prefs } from "./prefs.js";

mountAccount();
const ui = mountApp(CHEM);
mountExtras();
mountStats("chemistry");
mountHome(CHEM, ui);
prefs.sync();
