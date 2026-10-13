import type { Metadata } from 'next';
import './globals.css';
export const metadata:Metadata={title:'Foil Studio — Airfoil workspace',description:'A modern airfoil geometry and potential-flow analysis workspace.'};
export default function Layout({children}:{children:React.ReactNode}){return <html lang="en"><body>{children}</body></html>}
